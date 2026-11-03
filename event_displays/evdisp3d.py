"""Interactive 3D event viewer (three.js, schematic geometry) -> out/evdisp/bic_event_viewer.html
Fills evdisp3d_template.html with out/evdisp/events.json (written by evdisp.py), rounded to 0.1."""
import json
d = json.load(open('out/evdisp/events.json'))
r = lambda x: round(x, 1) if isinstance(x, float) else [r(v) for v in x] if isinstance(x, list) else {k: r(v) for k, v in x.items()} if isinstance(x, dict) else x
for e in d:
    e.pop('sfil5', None); e['astro']['purity'] = round(e['astro']['purity'], 3)
html = open('evdisp3d_template.html').read().replace('__DATA__', json.dumps(r(d), separators=(',', ':')))
open('out/evdisp/bic_event_viewer.html', 'w').write(html); print("wrote out/evdisp/bic_event_viewer.html")
