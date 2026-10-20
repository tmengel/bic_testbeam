"""Shared ROOT style for the run-3062 figures.

The sPHENIX plotting style (sPhenixStyle.C from the sPHENIX rootmacros) with no experiment label.
That style switches histogram titles off, so finish(canvas) writes each pad's description (the
title part of the first histogram or graph drawn in it) as a label in the top margin, outside the
frame: descriptions never cover data. Legends come from legend(): no box, no fill. Frames that hold
a legend are given headroom with headroom()."""
import ROOT
SPHENIX_STYLE = '/cvmfs/sphenix.sdcc.bnl.gov/alma9.2-gcc-14.2.0/release/release_new/new/rootmacros/sPhenixStyle.C'
TOP = 0.12                                   # top margin that holds the panel label

def setup(palette=None):
    ROOT.gROOT.SetBatch(True)
    ROOT.gErrorIgnoreLevel = ROOT.kWarning
    if not hasattr(ROOT, 'SetsPhenixStyle'):
        ROOT.gROOT.ProcessLine('std::cout.setstate(std::ios_base::failbit);')
        ROOT.gROOT.LoadMacro(SPHENIX_STYLE)
        ROOT.gROOT.ProcessLine('std::cout.clear();')
    ROOT.SetsPhenixStyle()
    s = ROOT.gStyle
    s.SetPadTopMargin(TOP); s.SetPadLeftMargin(0.17); s.SetPadBottomMargin(0.16); s.SetPadRightMargin(0.05)
    s.SetTitleOffset(1.3, "x"); s.SetTitleOffset(1.6, "y"); s.SetTitleOffset(1.3, "z")
    s.SetLegendFont(42); s.SetLegendFillColor(0); s.SetLegendBorderSize(0)
    s.SetPaintTextFormat(".2f")
    ROOT.TGaxis.SetMaxDigits(6)                # no x10^n exponent colliding with the panel label
    s.SetPalette(palette if palette is not None else ROOT.kBird)
    ROOT.gROOT.ForceStyle()

def colz(pad=None, right=0.17):
    """room on the right for a colour scale"""
    (pad or ROOT.gPad).SetRightMargin(right)

def legend(x1, y1, x2, y2, size=0.042, ncol=1):
    lg = ROOT.TLegend(x1, y1, x2, y2); lg.SetBorderSize(0); lg.SetFillStyle(0); lg.SetTextFont(42); lg.SetTextSize(size); lg.SetNColumns(ncol)
    return lg

def headroom(frame, ymin, ydata, frac=0.45, log=False):
    """set the frame's y range so the top `frac` of the plot is free for legends"""
    if log: frame.SetMaximum(ydata * 10 ** (frac / (1 - frac) * max(1.0, ROOT.TMath.Log10(ydata / max(ymin, 1e-12)))))
    else: frame.SetMaximum(ymin + (ydata - ymin) / (1 - frac))
    frame.SetMinimum(ymin)

def label(text, x=None, y=None, size=0.045, pad=None, align=11, color=ROOT.kBlack):
    """NDC text; default: top-left inside the top margin"""
    p = pad or ROOT.gPad
    t = ROOT.TLatex(); t.SetNDC(); t.SetTextFont(42); t.SetTextSize(size); t.SetTextAlign(align); t.SetTextColor(color)
    t.DrawLatex(p.GetLeftMargin() if x is None else x, 1 - p.GetTopMargin() + 0.025 if y is None else y, text)
    return t

def _title_of(pad):
    for o in pad.GetListOfPrimitives():
        if isinstance(o, (ROOT.TH1, ROOT.TGraph, ROOT.TMultiGraph)) and o.GetTitle():
            return o.GetTitle()
    return ""

def _fit(text, size, p):
    """shrink a label so it fits between the pad's left edge of the frame and its right border"""
    import re
    n = len(re.sub(r'#[a-zA-Z]+|[{}_^]', '', text))
    w, h = p.GetWw() * p.GetAbsWNDC(), p.GetWh() * p.GetAbsHNDC()
    return min(size, (1 - p.GetLeftMargin() - 0.02) / (0.5 * max(n, 1) * h / max(w, 1)))

def finish(canvas, size=0.045):
    """write every pad's description into its top margin (sPHENIX style has no titles)"""
    pads = [canvas.GetPad(i) for i in range(1, 100) if canvas.GetPad(i)] or [canvas]
    for p in pads:
        t = _title_of(p)
        if t:
            p.cd(); label(t, size=_fit(t, size, p), pad=p)
        p.Modified()
    canvas.Update()
