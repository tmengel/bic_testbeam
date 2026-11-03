"""bic_streams: read the run-3062 stream trees (KCU0, KCU1, CAEN, AstroPix) together, match triggers
between streams and work on waveforms (pedestal subtraction, pulse fits). Meant to be imported and
played with; `python3 bic_streams.py` runs a short demo.

    from bic_streams import Streams, pedestal, pulse_features, fit_pulse
    S = Streams()                      # files from bic_testbeam/data/*_3062.root (symlinks to tg01)
    S.index()                          # ts0 / gtt / counters of every trigger (cached in common/cache/)
    j1 = S.match_kcu1()                # KCU0 entry -> KCU1 entry (-1 if none), exact tick offset
    jc = S.match_caen()                # KCU0 entry -> CAEN entry, linear clock model + running median
    ev = S.event(123456)               # waveforms of one KCU0 trigger with its partners
    ped, wf = pedestal(ev['kcu0']['sfil3']['Nc3_adc'])

Clock constants (docs/method.md, docs/root_files.md):
  KCU tick 6.25 ns; ts(KCU1) = ts(KCU0) + KCU1_OFFSET (+0 or +4 ticks);
  KCU0 ticks = CAEN_RATIO x gtt + b(t), b found by a coarse scan then a running median (3 us wander);
  AstroPix on-chip ts = ((K + 4) // 16 - ASTRO_PHASE) mod 256 for a KCU0 trigger at tick K;
  readout stamp fpga_ts (25 ns) lags the trigger by 19 or 88 ms (+-5 ms).
"""
import os, sys, numpy as np, uproot
try:
    from scipy.optimize import curve_fit
except Exception:
    curve_fit = None

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, '..', 'data')
TICK_NS = 6.25
KCU1_OFFSET = -1143665946           # ts(KCU1) - ts(KCU0), run 3062
CAEN_RATIO = 1.365340679            # KCU ticks per CAEN GROUPTRIGGERTIME tick (8.533 ns)
ASTRO_PHASE = 189                   # ((K + 4) // 16 - ts) mod 256 == ASTRO_PHASE for true partners
ASTRO_LATENCY_MS = (19., 88.)       # trigger -> fpga_ts readout stamp
SFIL_S = {1: [2, 6, 10, 14], 2: [21, 25, 29, 33], 3: [63, 59, 52, 48], 4: [43, 37, 0, 4]}
SFIL_N = {1: [146, 150, 154, 158], 2: [165, 169, 173, 177], 3: [207, 203, 196, 192], 4: [187, 186, 144, 148]}
PADDLES = {'P3': 84, 'P4': 86, 'P5': 81, 'P6': 97}

TOA_WRAP_NS = 6.37                  # the LED macro's "magic number" (its configuration); only used in macro mode
def toa_ns(code):
    """nominal TOA code -> ns: fine[2:0] x 25 ps + mid[7:3] x 0.2 ns + coarse[9:8] x 6.25 ns (non-linear in practice)"""
    code = np.asarray(code, dtype=np.int64); return (code & 7) * 0.025 + ((code >> 3) & 31) * 0.2 + ((code >> 8) & 3) * 6.25
def toa_time(sample, code, macro=False, wrap=TOA_WRAP_NS):
    """time (ns from the start of sample 0) of the TOA word found in `sample`.
    Default (what run 3062 supports, method.md 4.9): 25 x sample + TOA_ns, no window rule.
    macro=True: the LED macro's rule, 25 x (sample + (0 if TOA_ns < wrap else 1)) - TOA_ns."""
    t = toa_ns(code); sample = np.asarray(sample)
    return (np.where(t < wrap, sample, sample + 1) * 25.0 - t) if macro else sample * 25.0 + t

# ---------------------------------------------------------------- waveform tools (KCU: 20 samples of 25 ns, 10 bit; CAEN: 256 of 1 ns, 12 bit, negative)
def pedestal(wf, n=4, axis=-1):
    """pedestal from the first n samples; returns (pedestal, pedestal-subtracted waveform)"""
    wf = np.asarray(wf, dtype=float); ped = wf[..., :n].mean(axis=axis, keepdims=True)
    return ped.squeeze(axis), wf - ped

def pulse_features(wf, n_ped=4, frac=0.2, negative=False, sat=1023):
    """pedestal, amplitude (max above pedestal), crossing time at `frac` of the amplitude (linear interpolation,
    in samples), integral, and whether any sample is saturated. negative=True for the paddles and the CAEN."""
    wf = np.asarray(wf, dtype=float); s = -wf if negative else wf
    ped, v = pedestal(s, n_ped); v = np.atleast_2d(v); ped = np.atleast_1d(ped)
    amp = v.max(axis=1); imax = v.argmax(axis=1); thr = frac * amp
    tcross = np.full(len(v), np.nan)
    for i in range(len(v)):
        if amp[i] <= 0: continue
        above = np.flatnonzero(v[i, :imax[i] + 1] >= thr[i])
        if len(above) == 0 or above[0] == 0: continue
        k = above[0]; y0, y1 = v[i, k - 1], v[i, k]
        tcross[i] = k - 1 + (thr[i] - y0) / (y1 - y0) if y1 != y0 else k
    saturated = (wf == (0 if negative else sat)).any(axis=-1) if wf.ndim > 1 else bool((wf == (0 if negative else sat)).any())
    return dict(ped=ped.squeeze(), amp=amp.squeeze(), tcross=tcross.squeeze(), integral=v.sum(axis=1).squeeze(), saturated=saturated)

def pulse_model(t, A, t0, tau_r, tau_f, p):
    """A x ((t - t0) / tau_r)^p x exp(-(t - t0 - p tau_r) / tau_f) for t > t0 (peak at t0 + p tau_r)"""
    x = np.clip((t - t0) / tau_r, 0, None)
    return A * x ** p * np.exp(-(t - t0 - p * tau_r) / tau_f) * (t > t0)

def fit_pulse(wf, dt=25., n_ped=4, sat=1023, negative=False):
    """fit one waveform with pulse_model, ignoring saturated samples; returns (params, pedestal, mask) or None.
    Time unit: ns (dt = sample spacing). Params: A, t0, tau_r, tau_f, p."""
    if curve_fit is None: raise RuntimeError('scipy is needed for fit_pulse')
    wf = np.asarray(wf, dtype=float); s = -wf if negative else wf
    ped, v = pedestal(s, n_ped); t = np.arange(len(v)) * dt
    mask = (wf != (0 if negative else sat))
    if mask.sum() < 6 or v.max() < 20: return None
    imax = v.argmax(); p0 = [v.max(), max(t[imax] - 2 * dt, 0), dt, 4 * dt, 2.]
    try:
        par, cov = curve_fit(pulse_model, t[mask], v[mask], p0=p0, bounds=([0, -dt, 1, 1, 0.5], [1e5, t[-1], 10 * dt, 50 * dt, 8]), maxfev=4000)
    except Exception:
        return None
    return dict(par=par, ped=ped, mask=mask, t=t, v=v, model=pulse_model(t, *par))

# ---------------------------------------------------------------- the streams
class Streams:
    def __init__(self, data=DATA, cache=os.path.join(HERE, 'cache')):
        self.path = {k: os.path.join(data, f'{k}_3062.root') for k in ('kcu0', 'kcu1', 'caen', 'astropix')}
        self.f = {k: uproot.open(p) for k, p in self.path.items() if os.path.exists(p)}
        self.cache = cache; os.makedirs(cache, exist_ok=True); self.idx = {}

    # ---- per-trigger indices (small arrays over the whole run, cached)
    def index(self, force=False):
        c = os.path.join(self.cache, 'index.npz')
        if os.path.exists(c) and not force:
            d = np.load(c); self.idx = {k: d[k] for k in d.files}; return self.idx
        h = self.f['kcu0']['header'].arrays(['ts0', 'complete', 'rt_first'], library='np')
        self.idx['k0_ts'] = h['ts0']; self.idx['k0_complete'] = h['complete']; self.idx['k0_rt'] = h['rt_first']
        if 'kcu1' in self.f:
            h = self.f['kcu1']['header'].arrays(['ts0', 'complete'], library='np'); self.idx['k1_ts'] = h['ts0']; self.idx['k1_complete'] = h['complete']
        if 'caen' in self.f:
            h = self.f['caen']['header'].arrays(['gtt', 'event_counter', 'rt'], library='np'); self.idx['c_gtt'] = h['gtt']; self.idx['c_evnr'] = h['event_counter']; self.idx['c_rt'] = h['rt']
        if 'astropix' in self.f:
            h = self.f['astropix']['hits'].arrays(['fpga_ts', 'ts', 'id', 'valid'], library='np')
            f = h['fpga_ts'].astype(np.int64); wraps = np.cumsum(np.diff(f, prepend=f[0]) < -2**31); self.idx['a_fpga'] = f + wraps * 2**32
            self.idx['a_ts'] = h['ts']; self.idx['a_chip'] = h['id']; self.idx['a_valid'] = h['valid']
        np.savez(c, **self.idx); return self.idx

    # ---- KCU0 <-> KCU1: exact tick offset
    def match_kcu1(self, tol=8, offset=KCU1_OFFSET):
        """KCU0 entry -> KCU1 entry within tol ticks of ts0 + offset, else -1"""
        k0, k1 = self.idx['k0_ts'], self.idx['k1_ts']; order = np.argsort(k1); ks = k1[order]
        pred = k0 + offset; j = np.searchsorted(ks, pred); j = np.clip(j, 1, len(ks) - 1)
        near = np.where(np.abs(ks[j] - pred) < np.abs(ks[j - 1] - pred), j, j - 1)
        ok = np.abs(ks[near] - pred) <= tol; out = np.where(ok, order[near], -1)
        print(f'KCU1 partners: {ok.sum()} of {len(k0)} KCU0 triggers ({100 * ok.mean():.1f}%)')
        return out

    def find_offset(self, a, b, guess, window, fine_bin=1., coarse_bin=1000., n=5000):
        """offset b - a that maximises coincidences: coarse histogram of pairwise differences within +-window of guess, then fine"""
        a = a[:n]; b = b[(b > a[0] + guess - window) & (b < a[-1] + guess + window)][:n * 4]
        d = (b[None, :] - a[:, None]).ravel(); d = d[np.abs(d - guess) < window]
        h, e = np.histogram(d, bins=int(2 * window / coarse_bin), range=(guess - window, guess + window)); c0 = 0.5 * (e[h.argmax()] + e[h.argmax() + 1])
        d2 = d[np.abs(d - c0) < 2 * coarse_bin]; h, e = np.histogram(d2, bins=int(4 * coarse_bin / fine_bin), range=(c0 - 2 * coarse_bin, c0 + 2 * coarse_bin))
        return 0.5 * (e[h.argmax()] + e[h.argmax() + 1]), h.max() / max(np.median(h), 1)

    # ---- KCU0 <-> CAEN: linear clock model with a running-median correction
    def match_caen(self, tol_ns=100., ratio=CAEN_RATIO, block=201):
        k0 = self.idx['k0_ts'].astype(np.float64); gtt = self.idx['c_gtt'].astype(np.float64) * ratio   # CAEN in KCU ticks, no offset yet
        c = os.path.join(self.cache, 'caen_offset.npz')
        if os.path.exists(c): d = np.load(c); b0, slope = float(d['b0']), float(d['slope'])
        else:
            # coarse guess from the host clocks (CAEN lags up to 2 s), then coincidence scans at the start and the end of
            # the run: offset and a residual slope (the fixed ratio is good to ~0.1 ppm, which is still 100 us over the run)
            rt_guess = np.median(self.idx['k0_rt'][:20000] - k0[:20000] * TICK_NS) - np.percentile(self.idx['c_rt'][:20000] - gtt[:20000] * TICK_NS, 1)
            b0, sig = self.find_offset(gtt, k0, -rt_guess / TICK_NS, window=3e9 / TICK_NS, fine_bin=4., coarse_bin=160000.)
            b1, sig1 = self.find_offset(gtt[-6000:], k0[-30000:], b0, window=2e6 / TICK_NS, fine_bin=4., coarse_bin=800.)
            slope = (b1 - b0) / (gtt[-3000] - gtt[3000])
            print(f'CAEN offset: {b0 * TICK_NS * 1e-9:.6f} s at the start (peak / background {sig:.0f}), {(b1 - b0) * TICK_NS * 1e-3:+.1f} us later at the end (peak {sig1:.0f}) -> residual slope {slope * 1e6:+.3f} ppm'); np.savez(c, b0=b0, slope=slope)
        gtt = gtt + slope * (gtt - gtt[3000])
        pred = gtt + b0
        def nearest(pred):
            j = np.clip(np.searchsorted(k0, pred), 1, len(k0) - 1); near = np.where(np.abs(k0[j] - pred) < np.abs(k0[j - 1] - pred), j, j - 1); return near, (k0[near] - pred) * TICK_NS
        # the clock wanders by a few us around the line: follow it with a running median of the residual, twice
        for tol in (20000., 2000.):
            near, res = nearest(pred); idx = np.flatnonzero(np.abs(res) < tol)
            if len(idx) < block: break
            rm = np.array([np.median(res[idx[max(0, k - block // 2):k + block // 2 + 1]]) for k in range(len(idx))])
            pred = pred + np.interp(np.arange(len(pred)), idx, rm) / TICK_NS
        near, res2 = nearest(pred)
        good = np.abs(res2) < tol_ns
        out = np.full(len(k0), -1, dtype=np.int64); out[near[good]] = np.flatnonzero(good)
        print(f'CAEN partners: {good.sum()} of {len(gtt)} CAEN triggers within {tol_ns:.0f} ns ({100 * good.mean():.1f}%), residual MAD {np.median(np.abs(res2[good] - np.median(res2[good]))):.1f} ns')
        return out

    # ---- AstroPix: phase test for a KCU0 trigger
    @staticmethod
    def astro_phase(k0_ticks, ts):
        """0 for a true partner: ((K + 4) // 16 - ts - ASTRO_PHASE) mod 256, centred to -128..127"""
        return (((np.asarray(k0_ticks, dtype=np.int64) + 4) // 16 - np.asarray(ts, dtype=np.int64) - ASTRO_PHASE + 128) % 256) - 128

    # ---- one trigger with its partners
    def event(self, i, j1=None, jc=None, trees=('sfil1', 'sfil2', 'sfil3', 'sfil4', 'paddles')):
        out = {'i': int(i), 'kcu0': {}}
        for t in ('header',) + tuple(trees): out['kcu0'][t] = {k: v[0] for k, v in self.f['kcu0'][t].arrays(entry_start=i, entry_stop=i + 1, library='np').items()}
        if j1 is not None and j1[i] >= 0: out['kcu1'] = {t: {k: v[0] for k, v in self.f['kcu1'][t].arrays(entry_start=j1[i], entry_stop=j1[i] + 1, library='np').items()} for t in ('header', 'sfil5')}
        if jc is not None and jc[i] >= 0: out['caen'] = {t: {k: v[0] for k, v in self.f['caen'][t].arrays(entry_start=jc[i], entry_stop=jc[i] + 1, library='np').items()} for t in ('header', 'bbcal', 'paddles')}
        return out

    def sfil_waveforms(self, i, layer):
        """(S[4][20], N[4][20]) ADC arrays of one SFIL layer for KCU0 trigger i (cells 1-4; -1 = missing line)"""
        t = self.f['kcu0'][f'sfil{layer}']; a = t.arrays(entry_start=i, entry_stop=i + 1, library='np')
        return np.array([a[f'Sc{c}_adc'][0] for c in range(1, 5)]), np.array([a[f'Nc{c}_adc'][0] for c in range(1, 5)])

# ---------------------------------------------------------------- demo
if __name__ == '__main__':
    S = Streams(); S.index()
    print({k: v.shape for k, v in S.idx.items()})
    j1 = S.match_kcu1(); jc = S.match_caen()
    both = np.flatnonzero((j1 >= 0) & (jc >= 0) & (S.idx['k0_complete'] == 1)); print('KCU0 triggers with KCU1 and CAEN partners and a complete waveform:', len(both))
    i = int(both[len(both) // 2]); ev = S.event(i, j1, jc)
    s3, n3 = S.sfil_waveforms(i, 3)
    fs, fn = pulse_features(s3), pulse_features(n3)
    print(f'trigger {i}: SFIL3 South amp {np.round(fs["amp"], 0)}, North amp {np.round(fn["amp"], 0)}, crossing (samples) S {np.round(fs["tcross"], 2)} N {np.round(fn["tcross"], 2)}')
    fit = fit_pulse(n3[2]); print('North cell 3 fit:', None if fit is None else dict(zip(['A', 't0', 'tau_r', 'tau_f', 'p'], np.round(fit['par'], 2))))
    if 'caen' in ev: print('CAEN partner event_counter', ev['caen']['header']['event_counter'], 'paddle P4 min', ev['caen']['paddles']['P4_adc'].min())
    if 'astropix' in S.f:
        K = S.idx['k0_ts'][i]; tsel = (S.idx['a_fpga'] * 4 > K - 1) & S.idx['a_valid']   # placeholder window: all later hits; refine with the latency offset
        print('AstroPix hits whose on-chip phase matches this trigger (whole run, no latency cut):', int((S.astro_phase(K, S.idx['a_ts'][tsel]) == 0).sum()), 'of', int(tsel.sum()))
    # a picture
    try:
        import ROOT; sys.path.insert(0, HERE); from plotstyle import setup, label
        setup(); c = ROOT.TCanvas('c', '', 1600, 500); c.Divide(3, 1); keep = []
        c.cd(1); fr = ROOT.gPad.DrawFrame(-0.5, -50, 19.5, 1100, ';sample (25 ns);ADC - pedestal'); keep.append(fr)
        for k, (w, col, nm) in enumerate(((s3[2], 2, 'S3 cell 3'), (n3[2], 4, 'N3 cell 3'))):
            _, v = pedestal(w); g = ROOT.TGraph(20, np.arange(20, dtype=float), v); g.SetLineColor(col); g.SetMarkerColor(col); g.SetMarkerStyle(20); g.Draw('LP'); keep.append(g); label(nm, x=0.2, y=0.85 - 0.06 * k, color=col)
        if fit is not None: g = ROOT.TGraph(20, fit['t'] / 25., fit['model']); g.SetLineColor(4); g.SetLineStyle(2); g.Draw('L'); keep.append(g)
        c.cd(2); fr = ROOT.gPad.DrawFrame(-0.5, -50, 19.5, 1100, ';sample (25 ns);ADC - pedestal'); keep.append(fr)
        for L, col in ((1, 1), (2, 2), (3, 4), (4, 6)):
            _, nn = S.sfil_waveforms(i, L); _, v = pedestal(nn[2]); g = ROOT.TGraph(20, np.arange(20, dtype=float), v); g.SetLineColor(col); g.SetMarkerColor(col); g.SetMarkerStyle(20); g.Draw('LP'); keep.append(g); label(f'SFIL{L} North cell 3', x=0.2, y=0.85 - 0.06 * (L - 1), color=col)
        c.cd(3)
        if 'caen' in ev:
            fr = ROOT.gPad.DrawFrame(0, -3000, 256, 200, ';sample (1 ns);CAEN ADC - pedestal'); keep.append(fr)
            for k, (nm, col) in enumerate((('P3', 1), ('P4', 2), ('P5', 4), ('P6', 6))):
                _, v = pedestal(ev['caen']['paddles'][f'{nm}_adc'], n=20); g = ROOT.TGraph(256, np.arange(256, dtype=float), v); g.SetLineColor(col); g.Draw('L'); keep.append(g); label(f'CAEN {nm}', x=0.2, y=0.85 - 0.06 * k, color=col)
        c.SaveAs(os.path.join(HERE, 'demo_event.png')); print('wrote common/demo_event.png')
    except Exception as e:
        print('no figure:', e)
