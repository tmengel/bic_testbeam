"""Where the KCU0 TOA word sits relative to the ADC 50 % crossing, and how its code relates to the crossing phase.
    python3 toa_explore.py [nmax] -> toa_explore.txt/.png"""
import numpy as np, uproot, sys, ROOT
sys.path.insert(0, '/sphenix/user/tmengel/eic/bic_testbeam/common'); from plotstyle import setup
setup(); ROOT.gROOT.SetBatch(True); out = open('toa_explore.txt', 'w')
def say(*a): print(*a, flush=True); print(*a, file=out)
NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 600000
f = uproot.open('/sphenix/user/tmengel/eic/bic_testbeam/data/kcu0_3062.root')
def load(tree, br, neg=False):
    a = f[tree].arrays([br + '_adc', br + '_toa'], entry_stop=NMAX, library='np'); W = a[br + '_adc'].astype(float); T = a[br + '_toa'].astype(int)
    good = W[:, 0] >= 0; ped = W[:, 1:4].mean(axis=1); y = W - ped[:, None]
    if neg: y = -y
    ipk = y[:, 4:16].argmax(axis=1) + 4; n = len(y); A = y[np.arange(n), ipk]
    half = 0.5 * A; k = ipk.copy()
    for d in range(1, 6):
        kk = np.clip(ipk - d, 0, 19); k = np.where(y[np.arange(n), kk] >= half, kk, k)
    k = np.clip(k, 1, 19); y0 = y[np.arange(n), k - 1]; y1 = y[np.arange(n), k]; frac = np.clip((half - y0) / np.maximum(y1 - y0, 1e-3), 0, 1)
    tc = 25. * (k - 1 + frac)   # crossing time; crossing lies inside sample interval [k-1, k]
    return y, T, good, A, ipk, k, tc
CH = [('SFIL2 S c3', 'sfil2', 'Sc3', False), ('SFIL2 N c3', 'sfil2', 'Nc3', False), ('P3', 'paddles', 'P3', True), ('P4', 'paddles', 'P4', True), ('P5', 'paddles', 'P5', True), ('P6', 'paddles', 'P6', True)]
c = ROOT.TCanvas('c', '', 2400, 1600); c.Divide(3, 2); keep = []
for i, (name, tree, br, neg) in enumerate(CH):
    y, T, good, A, ipk, k, tc = load(tree, br, neg)
    lo, hi = (200, 900) if not neg else (150, 1e9)
    sel = good & (A > lo) & (A < hi) & (ipk >= 5) & (ipk <= 14)
    n = sel.sum(); say('\n%s: %d pulses (A %g-%g)' % (name, n, lo, hi))
    ped_rate = ((T[good][:, 1:4] > 0).mean()); say('  word rate on pedestal samples 1-3: %.3f per sample' % ped_rate)
    for d in range(-2, 4):
        j = np.clip(k + d, 0, 19); w = T[np.arange(len(T)), j]; say('  word on crossing sample %+d: %.3f' % (d, (w[sel] > 0).mean()))
    # words at k+1 (expected): code vs phase within the sample before the word
    h = ROOT.TH2F('h%d' % i, '%s;t_{cross} - 25 (s_{word} - 1) [ns];TOA code' % name, 100, -25, 50, 128, 0, 1024)
    for d in (0, 1, 2):
        j = np.clip(k + d, 0, 19); w = T[np.arange(len(T)), j]; m = sel & (w > 0)
        ph = tc[m] - 25. * (j[m] - 1)
        for a_, b_ in zip(ph, w[m]): h.Fill(a_, b_)
    c.cd(i + 1); ROOT.gPad.SetLogz(); h.Draw('colz'); keep.append(h)
    # tightness: words at k+1, per code bin of 32: median phase and MAD
    j = np.clip(k + 1, 0, 19); w = T[np.arange(len(T)), j]; m = sel & (w > 0); ph = tc[m] - 25. * (j[m] - 1); cw = w[m]
    rows = []
    for b in range(32):
        mm = (cw // 32) == b
        if mm.sum() > 100: med = np.median(ph[mm]); rows.append('%4d:%5.1f/%4.1f(%d)' % (b * 32, med, 1.4826 * np.median(np.abs(ph[mm] - med)), mm.sum()))
    say('  words at k+1, code-bin: median phase / MAD-sigma (n):'); [say('   ' + '  '.join(rows[r:r + 4])) for r in range(0, len(rows), 4)]
c.SaveAs('toa_explore.png'); say('wrote toa_explore.png')
