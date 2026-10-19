# Sagnac effect: slide images

This folder holds 14 main slides and 3 backup slides. Each PNG is a full 1920 × 1080 (16:9) slide. Slides 02 and 06 also have an animation (`*_ANIMATION.gif`) that goes on top of the base slide. Slides 10 (random pixels) and 12 (drift) are optional and can be skipped if time is short.

Speaker notes for every slide, including the source document of each formula, are in **`presenter_notes.md`**.

All numbers come from the linear-calibration analysis in `../sagnac_update/` (`sagnac.json`, `sagnac.txt`, `sagnac_speeds.csv`, `methods.json`). The photographs, ring fits and stopwatch data come from the main folder. Only backup B1 uses the quadratic calibration.

## Putting them into Google Slides

1. **File → Page setup → Widescreen 16:9.** This is 10 × 5.625 in (25.4 × 14.29 cm).
2. For each PNG, choose **Insert → Image** and drag it to fill the slide. You can also set **Format options → Size & rotation**: width 10 in, then **Position**: X 0, Y 0.
3. For the two animations, insert the GIF on top of its base slide and type in this size and position under **Format options**:

| Slide | Animation | X | Y | Width | Height |
|---|---|---|---|---|---|
| 02 | `02_sagnac_principle_ANIMATION.gif` | 0.5 in (1.27 cm) | 1.0 in (2.54 cm) | 4.5 in (11.43 cm) | 4.0 in (10.16 cm) |
| 06 | `06_rotation_sweep_ANIMATION.gif` | 0.5 in (1.27 cm) | 1.0 in (2.54 cm) | 9.0 in (22.86 cm) | 4.0 in (10.16 cm) |

If your page has a different size, multiply all four numbers by (your slide width ÷ 10 in). Each base PNG already shows the animation's last frame in that rectangle, so the slide still reads correctly in a PDF export or a printout.

- Animation 02 loops every 7.3 s. Two pulses travel round a rotating loop, and the bar at the bottom marks when each returns to the splitter.
- Animation 06 loops every 20.4 s. It shows all 161 photographs in recording order, 0.1 s each, with the fitted dark rings drawn on. Δφ against ω fills in alongside, and the recording clock runs underneath. At the end it holds for 3.5 s on the fitted line and the prediction.

## Rebuilding

The scripts are in `source/`. From this folder, run:

```
python source/make_all.py
```

You need numpy, scipy, pandas, matplotlib and Pillow (the chaotic-circuit `.venv` has them). With ffmpeg on the PATH the GIFs come out smaller. Without it, Pillow writes them.
