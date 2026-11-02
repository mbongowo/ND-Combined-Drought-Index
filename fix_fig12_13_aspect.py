"""Restore North Dakota's true proportions in the land-cover map panels of Figures 12 and 13.

The map panels were drawn compressed east-west (width/height 1.22-1.56). Each panel is
resampled to the 1.645 width/height ratio of the state in NAD83 / UTM zone 14N (the same
projection as Figure 6) and re-centred in its slot. The donut charts are untouched.

usage: python fix_fig12_13_aspect.py IN_PNG OUT_PNG
"""
import sys
import numpy as np
from PIL import Image
from scipy import ndimage

ND_RATIO = 1.645
img = Image.open(sys.argv[1]).convert("RGB")
a = np.asarray(img).astype(int)
sat = (a.max(2) - a.min(2)) > 40
sat[:, int(a.shape[1] * 0.38):] = False
lab, _ = ndimage.label(ndimage.binary_closing(sat, iterations=3))
out = img.copy()
for sl in ndimage.find_objects(lab):
    y0, y1, x0, x1 = sl[0].start, sl[0].stop, sl[1].start, sl[1].stop
    if y1 - y0 < 100 or x1 - x0 < 150:
        continue
    m = 3
    y0, y1, x0, x1 = max(y0 - m, 0), y1 + m, max(x0 - m, 0), x1 + m
    panel = img.crop((x0, y0, x1, y1))
    w = x1 - x0; h = round(w / ND_RATIO)
    panel = panel.resize((w, h), Image.LANCZOS)
    out.paste((255, 255, 255), (x0, y0, x1, y1))
    out.paste(panel, (x0, y0 + ((y1 - y0) - h) // 2))
    print(f"panel {x0},{y0}: {w}x{y1-y0} -> {w}x{h}")
out.save(sys.argv[2], dpi=(300, 300))
