"""Rebuild every Sagnac slide image and both animations.

    python make_all.py

Needs numpy, scipy, pandas, matplotlib and Pillow (ffmpeg optional, for smaller GIFs). Reads the analysis
outputs in the 'Sagnac effect' folder (sagnac_update/ for the linear-calibration results), the photographs,
and random_points_scan.json (made by random_points_scan.py, about an hour on two cores).
"""
import slides
import anim_concept
import anim_measurements

for name in sorted(slides.REGISTRY, key=lambda n: (n[0] == 'b', n)):     # main slides, then backups
    slides.REGISTRY[name]()
anim_concept.main()
anim_measurements.main()
