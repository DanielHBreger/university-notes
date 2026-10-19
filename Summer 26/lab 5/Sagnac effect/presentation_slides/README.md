# Sagnac effect: slide images

This folder holds 13 main slides and 4 backup slides (B1–B4), and `overview.png` with all of them at a glance. Each PNG is a full 1920 × 1080 (16:9) slide. Slides 02 and 05 also have an animation (`*_ANIMATION.gif`) that goes on top of the base slide. Slide 10 (random pixels) is optional and can be skipped if time is short.

Speaker notes for every slide, including the source document of each formula, are in **`presenter_notes.md`** (short, to keep open while presenting). **`presenter_guide.md`** is the full guide for whoever presents: background physics and statistics, every slide explained element by element with a script, where each number comes from, a question bank and a numbers cheat sheet.

Every number on the slides and in the documents is computed by `source/collect_numbers.py`, which writes `source/numbers.json`. The results documents for the whole experiment are `../RESULTS.md` and `../SUMMARY.md`.

All numbers come from the linear-calibration analysis in `../sagnac_update/` (`sagnac.json`, `sagnac.txt`, `sagnac_speeds.csv`, `methods.json`). The photographs, ring fits and stopwatch data come from the main folder. Slide 10 uses `source/random_points_scan.json`, the random-pixel method of `methods.py` repeated for 1000 draws. The deck refits every phase-reading method with the same equal-weight linear regression and adds the calibration-model term to the uncertainty, so its total for the slope (± 0.07) is slightly larger than the report's (± 0.06). The quadratic calibration is used only for comparison: backup B1 shows it, and the 1.7 % difference between the linear and the quadratic slope is part of the uncertainty on slides 9, 11 and 12.

## Putting them into Google Slides

1. **File → Page setup → Widescreen 16:9.** This is 10 × 5.625 in (25.4 × 14.29 cm).
2. For each PNG, choose **Insert → Image** and drag it to fill the slide. You can also set **Format options → Size & rotation**: width 10 in, then **Position**: X 0, Y 0.
3. Slides 02 and 05 show a dashed gray box with the name of a GIF and four numbers. Insert that GIF on top of the slide and type the box's size and position under **Format options**. The GIF covers the box exactly.

| Slide | Animation | X | Y | Width | Height |
|---|---|---|---|---|---|
| 02 | `02_sagnac_principle_ANIMATION.gif` | 0.5 in (1.27 cm) | 1.0 in (2.54 cm) | 4.5 in (11.43 cm) | 4.0 in (10.16 cm) |
| 05 | `05_measurements_ANIMATION.gif` | 0.5 in (1.27 cm) | 1.0 in (2.54 cm) | 9.0 in (22.86 cm) | 4.0 in (10.16 cm) |

If your page has a different size, multiply all four numbers by (your slide width ÷ 10 in).

- Animation 02 loops every 7.3 s. Two pulses travel round a rotating loop, and the bar at the bottom marks when each returns to the splitter.
- Animation 05 loops every 19.4 s. It shows all 161 photographs in recording order, 0.1 s each, with the first dark ring at rest as a dashed reference circle. ω, the motor setting and the camera clock run along the top, and the rate of rotation during the run fills in on the right.

## Rebuilding

The scripts are in `source/`. From this folder, run:

```
python source/make_all.py
```

You need numpy, scipy, pandas, matplotlib and Pillow (the chaotic-circuit `.venv` has them). With ffmpeg on the PATH the GIFs come out smaller. Without it, Pillow writes them.

`source/random_points_scan.json` is already included. To remake it, run `python source/random_points_scan.py` (about an hour on two cores). It reuses the functions of `sagnac_update/methods.py` unchanged.
