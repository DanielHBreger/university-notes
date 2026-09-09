"""Index every final scientific PNG and verify its file integrity."""
from pathlib import Path
from PIL import Image,ImageStat,ImageOps,ImageDraw
import json,hashlib,html
from urllib.parse import quote
ROOT=Path(__file__).resolve().parents[1]
paths=sorted([*ROOT.glob('*.png'),*(ROOT/'chua').rglob('*.png'),
              *(ROOT/'verification').glob('M*.png'),*(ROOT/'verification').glob('trace*_M5.png')])
rows=[]
for p in paths:
    with Image.open(p) as im:
        im.verify()
    with Image.open(p) as im:
        assert min(im.size)>100 and max(ImageStat.Stat(im.convert('RGB')).stddev)>5,p
        rows.append(dict(path=p.relative_to(ROOT).as_posix(),size=list(im.size),
                         sha256=hashlib.sha256(p.read_bytes()).hexdigest(),integrity='passed'))
original=json.loads((ROOT/'verification/original_figures.json').read_text())
assert {r['path'] for r in original} <= {r['path'] for r in rows}
(ROOT/'verification/figure_manifest.json').write_text(json.dumps(rows,indent=2))
cards=[]
for row in rows:
    p=row['path'];url=quote('../'+p,safe='/');category='Requirements' if p.startswith('verification/') else 'Simulation' if '/simulated_' in p else 'Lyapunov' if 'lyapunov' in p else 'Return maps' if 'lorenz' in p else 'Bifurcation' if 'bifurcation' in p else 'Current–voltage'
    cards.append(f'<article data-category="{category}"><a href="{url}" target="_blank"><img loading="lazy" src="{url}" alt="{html.escape(p)}"></a><p>{html.escape(p)}</p><small>{category} · {row["size"][0]} × {row["size"][1]}</small></article>')
page='''<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Chua lab · verified figures</title>
<style>body{font:16px system-ui,sans-serif;background:#f4f6f8;color:#142536;margin:0;padding:32px;max-width:1600px;margin:auto}h1{font-size:30px;margin:0 0 12px}p{line-height:1.5}select{font:inherit;padding:10px;border:1px solid #9aa7b2;border-radius:6px;margin:12px 0 25px}main{display:grid;grid-template-columns:repeat(auto-fit,minmax(400px,1fr));gap:20px}article{background:white;padding:15px;border:1px solid #d6dfe6;border-radius:8px}img{width:100%;height:300px;object-fit:contain}article p{overflow-wrap:anywhere;margin:8px 0;font-size:14px}small{color:#586776}a{color:#126b9c}@media(max-width:600px){body{padding:16px}main{grid-template-columns:1fr}}</style>
<h1>Chua lab · verified figures</h1><p>All original figures were regenerated. Click a figure to inspect its full resolution.<br>Raw measurements are preserved. Offset removal is a model assumption; the record labelled “double-scroll end” is nearly periodic.</p>
<p><a href="../VERIFICATION_REPORT.md">Verification report</a> · <a href="M4_M5_results.csv">Three-record numerical results</a> · <a href="figure_manifest.json">Figure manifest</a></p>
<label for="filter">Show </label><select id="filter"><option>All figures</option><option>Requirements</option><option>Current–voltage</option><option>Return maps</option><option>Bifurcation</option><option>Lyapunov</option><option>Simulation</option></select><main>'''+''.join(cards)+'''</main><script>document.querySelector('select').addEventListener('change',e=>{document.querySelectorAll('article').forEach(a=>a.hidden=e.target.value!=='All figures'&&a.dataset.category!==e.target.value)})</script></html>'''
(ROOT/'verification/figures.html').write_text(page,encoding='utf-8')
# Review sheets retain sufficient pixels for layout and axis checks.
preview=ROOT/'verification/after_previews';preview.mkdir(exist_ok=True)
selected=[ROOT/r['path'] for r in rows if 'lyapunov' in r['path'] or 'simulated_' in r['path'] or r['path'].startswith('verification/')]
for start in range(0,len(selected),4):
    sheet=Image.new('RGB',(2000,1500),'#eeeeee');draw=ImageDraw.Draw(sheet)
    for i,p in enumerate(selected[start:start+4]):
        x=(i%2)*1000;y=(i//2)*750
        with Image.open(p) as im:
            fitted=ImageOps.contain(im.convert('RGB'),(980,700))
            sheet.paste(fitted,(x+(1000-fitted.width)//2,y+30))
        draw.text((x+10,y+8),p.relative_to(ROOT).as_posix(),fill='black')
    sheet.save(preview/f'final_{start//4+1:02d}.jpg',quality=92)
print(f'{len(rows)} figures checked; all {len(original)} originals included.')
