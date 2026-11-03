"""Off-beam cell response relative to the beam cell from the template-fit trees (any run or part), the fitted-amplitude version
of cell_equalisation.py (which used raw peaks of the first 800k triggers of 3062).

    python3 cell_ratios_fits.py <fit file> [tag]   -> cell_ratios_fits_<tag>.txt

Per layer and paddle, over triggers whose beam cell (cell 3) has a fitted pulse with 100 < raw peak < 900 ADC (unsaturated):
the mean fitted amplitude of every cell (a cell without a fitted pulse counts 0: the pulse finder needs > 60 ADC, so the
off-beam means are lower limits of the light share; the fraction with a pulse is given too), divided by the mean of the beam
cell. South and North separately; k = 8 layer + 4 end + cell.
"""
import sys, numpy as np, uproot
F = sys.argv[1]; TAG = sys.argv[2] if len(sys.argv) > 2 else 'x'
t = uproot.open(F)['fits']; n = t.num_entries
S = np.zeros((5, 32)); S2 = np.zeros((5, 32)); NP = np.zeros((5, 32)); NT = np.zeros((5, 4)); OKB = np.zeros((5, 4))
for arr in t.iterate(['pad', 'ok', 'amp', 'afit'], step_size=300000, library='np'):
    pad, ok, amp, afit = arr['pad'], arr['ok'], arr['amp'], arr['afit']
    for L in range(4):
        kb = 8 * L + 2; sel = ok[:, kb] & (amp[:, kb] > 100) & (amp[:, kb] < 900)
        for ip, p in enumerate((0, 3, 4, 5, 6)):
            m = sel if p == 0 else sel & (pad == p)
            if not m.any(): continue
            NT[ip, L] += m.sum()
            a = np.where(ok[m], afit[m], 0.)
            S[ip, 8 * L:8 * L + 8] += a[:, 8 * L:8 * L + 8].sum(axis=0); S2[ip, 8 * L:8 * L + 8] += (a[:, 8 * L:8 * L + 8] ** 2).sum(axis=0); NP[ip, 8 * L:8 * L + 8] += ok[m][:, 8 * L:8 * L + 8].sum(axis=0)
out = open('cell_ratios_fits_%s.txt' % TAG, 'w')
def say(s): print(s); out.write(s + '\n')
say('%s: %d fit entries' % (F, n))
say('per layer and paddle: triggers with an unsaturated beam-cell pulse (100-900 ADC); then per cell <A_fit>/<A_beam> (fraction with a pulse)')
for ip, p in enumerate((0, 3, 4, 5, 6)):
    say('\npaddle %s' % ('any single' if p == 0 else 'P%d' % p))
    for L in range(4):
        nt = NT[ip, L]
        if nt < 100: continue
        line = '  SFIL%d  n = %8d  beam S %6.1f ADC |' % (L + 1, nt, S[ip, 8 * L + 2] / nt)
        for e, en in ((0, 'S'), (1, 'N')):
            kb = 8 * L + 4 * e + 2; mb = S[ip, kb] / nt
            cells = []
            for c in range(4):
                k = 8 * L + 4 * e + c; cells.append('%s%d %.3f (%.2f)' % (en, c + 1, S[ip, k] / nt / mb if mb > 0 else float('nan'), NP[ip, k] / nt))
            line += ' ' + '  '.join(cells) + ' |'
        say(line)
out.close()
