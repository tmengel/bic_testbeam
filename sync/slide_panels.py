"""Reduced figures for the slides, cut from the full figures (same pixels, no re-plotting).

Each full figure is a ROOT canvas divided into equal pads (columns x rows); a panel is pad (col, row).
  out/caen_kcu0_slide.png       caen_kcu0.png       3x1: panels (1,0) residual vs control, (2,0) zoom
  out/astro_summary_slide.png   astro_summary.png   3x2: (0,0) 10 MHz phase lock, (1,0) readout latency,
                                                         (1,1) matches per chip
  out/astro_geo_kcu0_slide.png  astro_geo_kcu0.png  3x2: (0,0) all clusters, (2,0) paddle hit vs x,
                                                         (1,1) SFIL2 vertical profile
Run after the full figures (reproduce.sh plots/diag)."""
from PIL import Image

def cut(src, grid, panels, dst):
    im = Image.open(src).convert('RGB'); W, H = im.size; nc, nr = grid; pw, ph = W / nc, H / nr
    tiles = [im.crop((int(c * pw), int(r * ph), int((c + 1) * pw), int((r + 1) * ph))) for c, r in panels]
    out = Image.new('RGB', (sum(t.size[0] for t in tiles), max(t.size[1] for t in tiles)), 'white'); x = 0
    for t in tiles: out.paste(t, (x, 0)); x += t.size[0]
    out.save(dst); print(dst, out.size)

cut('out/caen_kcu0.png', (3, 1), [(1, 0), (2, 0)], 'out/caen_kcu0_slide.png')
cut('out/astro_summary.png', (3, 2), [(0, 0), (1, 0), (1, 1)], 'out/astro_summary_slide.png')
cut('out/astro_geo_kcu0.png', (3, 2), [(0, 0), (2, 0), (1, 1)], 'out/astro_geo_kcu0_slide.png')
