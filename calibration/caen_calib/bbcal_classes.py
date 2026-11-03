"""Baby BCAL South after the DRS4 cell-pedestal correction: channel shares in showers, SFIL-quiet triggers (what the
Baby BCAL sees when the SFILs see nothing), MIP-track search by column, and the cosmic coincidence budget.
Join: caen_calib_seg*.npz (global CAEN ordinal = seg offset + idx) -> triggers_3062.npz i_caen -> events_3062 rows.
    python3 bbcal_classes.py  -> bbcal_classes.txt, bbcal_classes.png
"""
import numpy as np, uproot, sys, ROOT
sys.path.insert(0, '/sphenix/user/tmengel/eic/bic_testbeam/common'); from plotstyle import setup, legend, label
setup(); out = open('bbcal_classes.txt', 'w')
def say(*a): print(*a, flush=True); print(*a, file=out)
import os
NSEG = [687701, 694330, 125979]; OFF = np.cumsum([0] + NSEG); NC = OFF[-1]
peak = np.full((NC, 16), np.nan, np.float32); tpk = np.zeros((NC, 16), np.int16); qn = np.full((NC, 16), np.nan, np.float32); q140 = np.full((NC, 16), np.nan, np.float32); rms0 = np.full((NC, 16), np.nan, np.float32)
have = np.zeros(NC, bool)
for s in (0, 1, 2):
    fn = '/sphenix/tg/tg01/jets/tmengel/eic/caen/calib/caen_calib_seg%d.npz' % s
    if not os.path.exists(fn): say('segment %d not available yet' % s); continue
    d = np.load(fn); sl = slice(OFF[s], OFF[s + 1]); assert len(d['peak']) == NSEG[s]
    peak[sl], tpk[sl], qn[sl], q140[sl], rms0[sl] = d['peak'], d['tpk'], d['qn'], d['q140'], d['rms0']; have[sl] = True
say('CAEN events with corrected data %d of %d; corrected baseline rms per channel (mean): %s' % (have.sum(), NC, np.nanmean(rms0, axis=0).round(1)))
TT = np.load('/sphenix/user/tmengel/eic/bic_testbeam/sync/out/triggers_3062.npz'); ic = TT['i_caen'].copy(); ic[(ic >= 0) & ~have[np.maximum(ic, 0)]] = -1
ev = uproot.open('/sphenix/user/tmengel/eic/bic_testbeam/data/events_3062.root')['events'].arrays(['has_kcu0', 'has_caen', 'kcu0_frag', 'pad', 'sfil_s', 'caen_q'], library='np')
N = len(ev['pad']); assert N == len(ic), (N, len(ic))
hc = ic >= 0; say('event rows %d, with CAEN %d' % (N, hc.sum()))
# sanity: the old caen_q (20-159, raw pedestal) vs new q140 for the same rows
old = ev['caen_q'][hc][:, :12].astype(float); new = q140[ic[hc]][:, :12]
say('old caen_q vs new q140 (same window, cell-corrected): per-channel correlation %s; mean old %s; mean new %s' % (np.round([np.corrcoef(old[:, i], new[:, i])[0, 1] for i in range(12)], 3), old.mean(axis=0).round(0), new.mean(axis=0).round(0)))
labels = ['A0', 'A1', 'A3', 'A4', 'A2', 'A5', 'A6', 'A7', 'A8', 'A9', 'A14', 'A11']   # run-log sheet 'Changed on Aug 10': S A0 0, S A1 1, S A3 2, S A4 3, (ch 4 unlabeled: A2 by its pattern), S A5 5, 6-9 -> A6-A9, 14 -> 10, 11 -> 11
col = {l: int(l[1:]) % 4 for l in labels}; narr = [int(l[1:]) // 4 + 1 for l in labels]
L = np.clip(ev['sfil_s'], 0, None).astype(float).sum(axis=2); S = L.sum(axis=1)
base = (ev['has_kcu0'] == 1) & hc & (ev['kcu0_frag'] == 0) & np.isin(ev['pad'], [1, 2, 4, 8])
Q = np.zeros((N, 16), np.float32); Q[hc] = qn[ic[hc]]; PK = np.zeros((N, 16), np.float32); PK[hc] = peak[ic[hc]]; TP = np.zeros((N, 16), np.int16); TP[hc] = tpk[ic[hc]]
Qb = np.clip(Q[:, :12], 0, None); PKb = PK[:, :12]
shower = base & (S > 600)
v = Qb[shower].mean(axis=0); say('\nshower shares (qn, corrected), all paddles: ' + ' '.join('%s %.3f' % (l, x) for l, x in zip(labels, v / v.sum())))
PB = {4: 2, 3: 1, 6: 8, 5: 4}
for p in (4, 3, 6, 5):
    v = Qb[shower & (ev['pad'] == PB[p])].mean(axis=0); say('  P%d: ' % p + ' '.join('%.3f' % x for x in v / v.sum()))
c = np.corrcoef(Qb[shower].T); say('shower correlation matrix (corrected):'); [say('  %-4s ' % labels[i] + ' '.join('%5.2f' % c[i, j] for j in range(12))) for i in range(12)]
# ---- noise threshold: 5 sigma on the peak (sigma ~ rms0), and the on-count distributions
thr = 5 * np.nanmean(rms0, axis=0)[:12]; on = PKb > thr[None, :]; non = on.sum(axis=1)
say('\npeak threshold per channel (5 x baseline rms): %s' % thr.round(0))
quiet = base & (S < 100); mid = base & (S >= 100) & (S < 300)
for name, m in (('showers', shower), ('SFIL quiet (S<100)', quiet), ('SFIL 100-300', mid)):
    say('%s: %d triggers; channels above threshold: %s' % (name, m.sum(), np.bincount(non[m], minlength=13)))
# quiet triggers with something in the Baby BCAL: pattern concentration and column structure
say('fraction of triggers with each channel above threshold, showers / quiet / 100-300:'); [say('  %-4s %.3f / %.3f / %.3f' % (labels[i], on[shower, i].mean(), on[quiet, i].mean(), on[mid, i].mean())) for i in range(12)]
qb = quiet & (non >= 1)
conc = Qb.max(axis=1) / np.maximum(Qb.sum(axis=1), 1)
say('\nSFIL-quiet triggers with >= 1 channel on: %d; of these max-channel share > 0.5: %d, > 0.8: %d' % (qb.sum(), (qb & (conc > 0.5)).sum(), (qb & (conc > 0.8)).sum()))
big = quiet & (Qb.sum(axis=1) > 20000)
say('SFIL-quiet with Baby BCAL sum > 20000 (direct hits?): %d; mean share pattern: %s' % (big.sum(), ' '.join('%s %.3f' % (l, x) for l, x in zip(labels, Qb[big].mean(axis=0) / Qb[big].mean(axis=0).sum()))))
say('   their peak-sample mean per channel (showers: %s): %s' % (np.round(TP[shower][:, :12].mean(axis=0), 0), np.round(TP[big][:, :12].mean(axis=0), 0)))
# ---- MIP-track search: triggers (any SFIL class) where exactly the channels of one column are on and nothing else
cols = {cc: [k for k in range(12) if col[labels[k]] == cc] for cc in range(4)}
say('\ncolumn-track search (all channels of one column on, all others off), by SFIL class:')
for cc in range(4):
    ks = cols[cc]; others = [k for k in range(12) if k not in ks]; track = on[:, ks].all(axis=1) & ~on[:, others].any(axis=1)
    for name, m in (('quiet', quiet), ('100-300', mid), ('showers', shower)):
        s = track & m
        if s.sum() == 0: say('  column %d (%s) %-8s: 0' % (cc, ','.join(labels[k] for k in ks), name)); continue
        say('  column %d (%s) %-8s: %5d triggers; median charge per channel %s; per array %s' % (cc, ','.join(labels[k] for k in ks), name, s.sum(), np.round(np.median(Qb[s][:, ks], axis=0), 0), np.round(np.median(Qb[s][:, ks], axis=0) / np.array([narr[k] for k in ks]), 0)))
# ---- cosmic budget
say('\ncosmic coincidence budget: ~1 muon / cm^2 / min through the 58 x 22 cm top face = ~21 Hz; CAEN window 256 ns; CAEN triggers %d -> expected coincidences ~ %.0f' % (len(peak), len(peak) * 21 * 256e-9))
# ---- figure: on-count distributions, concentration, column-track charge
c = ROOT.TCanvas('c', '', 2400, 800); c.Divide(3, 1); keep = []
c.cd(1); fr = ROOT.gPad.DrawFrame(-0.5, 0.5, 12.5, 1e6, ';Baby BCAL channels above 5#sigma;triggers'); ROOT.gPad.SetLogy(); keep.append(fr)
for (name, m), colr in zip((('showers', shower), ('SFIL quiet', quiet), ('SFIL 100-300', mid)), (ROOT.kBlack, ROOT.kRed + 1, ROOT.kBlue + 1)):
    h = ROOT.TH1F('h%s' % name[:4], '', 13, -0.5, 12.5); [h.Fill(x) for x in non[m]]; h.SetLineColor(colr); h.SetLineWidth(2); h.Draw('hist same'); keep.append(h)
lg = legend(0.45, 0.7, 0.88, 0.88); [lg.AddEntry(keep[-3 + i], nm, 'l') for i, nm in enumerate(('showers', 'SFIL quiet', 'SFIL 100-300'))]; lg.Draw(); keep.append(lg)
c.cd(2); fr = ROOT.gPad.DrawFrame(0, 0.5, 1, 1e5, ';max-channel share of the Baby BCAL sum;SFIL-quiet triggers'); ROOT.gPad.SetLogy(); keep.append(fr)
h = ROOT.TH1F('hconc', '', 50, 0, 1); [h.Fill(x) for x in conc[qb]]; h.SetLineWidth(2); h.Draw('hist same'); keep.append(h)
c.cd(3); fr = ROOT.gPad.DrawFrame(0, 0, 60000, 1.0, ';Baby BCAL sum (qn) [counts x ns];fraction'); keep.append(fr)
for (name, m), colr in zip((('showers', shower), ('SFIL quiet', quiet)), (ROOT.kBlack, ROOT.kRed + 1)):
    h = ROOT.TH1F('hs%s' % name[:4], '', 60, 0, 60000); [h.Fill(x) for x in Qb[m].sum(axis=1)]; h.Scale(1. / max(1, h.GetEntries())); h.SetLineColor(colr); h.SetLineWidth(2); h.Draw('hist same'); keep.append(h)
fr.SetMaximum(max(k.GetMaximum() for k in keep[-2:]) * 1.3)
c.SaveAs('bbcal_classes.png'); say('wrote bbcal_classes.png')
