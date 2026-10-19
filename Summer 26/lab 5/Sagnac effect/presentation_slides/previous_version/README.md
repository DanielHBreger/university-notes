# Sagnac effect: slide images

There are 12 main slides and 4 backup slides. Every PNG is a full 1920 × 1080 (16:9) slide. Slides 02 and 06 also have an animation (`*_ANIMATION.gif`) that goes on top of the base slide. The talk runs about 10 minutes.

Every number comes from the current linear-calibration analysis in `../sagnac_update/` (`sagnac.json`, `sagnac_speeds.csv`, `methods.json`). The photographs, ring fits and stopwatch data come from the main folder. Backup B1 is the only slide that uses the quadratic calibration.

## Putting them into Google Slides

1. **File → Page setup → Widescreen 16:9.** This is 10 × 5.625 in (25.4 × 14.29 cm).
2. For each PNG, choose **Insert → Image**, then drag the image to fill the slide. You can also use **Format options → Size & rotation**: width 10 in, then **Position**: X 0, Y 0.
3. For the two animations, insert the GIF on top of its base slide and type in this size and position under **Format options**:

| Slide | Animation | X | Y | Width | Height |
|---|---|---|---|---|---|
| 02 | `02_sagnac_principle_ANIMATION.gif` | 0.5 in (1.27 cm) | 1.0 in (2.54 cm) | 4.5 in (11.43 cm) | 4.0 in (10.16 cm) |
| 06 | `06_rotation_sweep_ANIMATION.gif` | 0.5 in (1.27 cm) | 1.0 in (2.54 cm) | 9.0 in (22.86 cm) | 4.0 in (10.16 cm) |

If your page has a different size, multiply all four numbers by (your slide width ÷ 10 in). Each base PNG already shows the animation's last frame in that rectangle, so the slide still reads correctly in a PDF export or a printout.

- Animation 02 loops every 7.3 s. Two pulses go round a rotating loop. The one moving with the rotation has to chase the moving splitter, so it arrives Δt later.
- Animation 06 loops every 20.4 s. It shows all 161 photographs in recording order, 0.1 s each, with the fitted dark rings drawn on. Δφ against ω fills in alongside, and a bar shows the recording clock. At the end it holds for 3.5 s on the fitted line and the prediction.

## Speaker notes (about 10 minutes)

**01 Title (0:20).** We measured the Sagnac effect: light sent both ways round a rotating fibre loop. These are our photographs. As the table turns faster, the rings move inward and the centre goes dark.

**02 Rotation delays the beam that travels with it (1:00). Let the animation run once.** Two beams leave the splitter together, and the splitter moves with the table while they travel. The co-rotating beam has to chase it and arrives later. The counter-rotating beam arrives early. For N turns of area A, the delay is Δt = 4NAω/c². Multiplying by the light's angular frequency gives Δφ = 8πNAω/(λc). At our top speed the delay is only 1.07 fs. That is half an optical period, so we expect about half a fringe. This is how fibre-optic gyroscopes work.

**03 Both beams travel the same loop (0:50).** The beam splitter sends one beam into the fibre at L₁ and the other at L₂. They go through the same 13 turns in opposite directions and recombine on the camera. Everything turns with the table. We measure two things independently: ω with a stopwatch, and Δφ from photographs. With N = 13 ± 1 and R = 57.0 cm, the prediction is a_th = 1.74 ± 0.14 rad per rad/s. The ±1 turn is almost all of that uncertainty.

**04 Rings, evenly spaced in r² (0:50).** The two beams leave the fibre with different wavefront curvature, so we see rings rather than straight fringes. Two spherical waves give rings that are evenly spaced in r², not in r. Here the spacing is 241 then 132 px in r, but 1.08 and 1.08 × 10⁵ px² in r². We therefore write the phase as 2πr²/P + φc. Rotation adds the same phase across the whole beam, which shifts φc. A tilt only moves the centre and a refocus only changes P, so φc is the right quantity to track.

**05 Calibration (0:45).** The motor setting S is not an angular velocity, so we timed 5–8 turns at nine settings. The manual's straight line gives k = 0.02039 rad/s per unit of S. The stopwatch is good to 0.1 %, but the table's speed itself varies by 0.48 %. The residuals show the table running 2.3 % below the line at the top speed. We keep the requested linear conversion; backup B1 shows the quadratic alternative.

**06 Faster rotation pulls the rings inward (0:45). Let it play once.** These are all 161 photographs in the order we took them. Watch the inner dark ring close in on the centre while the phase on the right climbs. Each small dot is one photo and each large dot is the mean at one setting. By the end, the points lie on a straight line.

**07 How a photo becomes a phase (0:50).** We find the dark rings and fit concentric circles. Then we plot R² against ring number, which gives a straight line with slope P and intercept η. The phase at the centre is φc = π − 2πη. From rest to full speed the intercept drops by 0.52 fringe, which is 3.3 rad. Only the ratio R²/P enters, so we need no pixel calibration. P is fitted for every photo, so a change in ring spacing does not bias the phase.

**08 The phase is linear in ω (1:00).** The slope is 1.84 ± 0.02 rad per rad/s, against a prediction of 1.74 ± 0.14, so they are 0.7σ apart. A quadratic term is consistent with zero. Between settings the phase scatters by 0.04 rad, more than the photo-to-photo noise explains. We add that scatter to the error bars as σx; without it, χ² = 57 for 8 degrees of freedom.

**09 Three ways of reading the phase agree to ±3 % (0:50).** Could the slope come from how we read the phase? We also used a cosine fit to each profile and template matching against the rest profile. The three slopes span 1.77–1.88. Half that range, 0.055, is our method uncertainty, so a = 1.84 ± 0.06. A fourth method, ten random pixels, gives a slope that depends on which pixels are drawn, so we leave it out (backup B3).

**10 The speed of light (0:45).** The same equation, solved for c, gives (2.8 ± 0.2) × 10⁸ m/s, 0.7σ from the defined value. The turn count supplies 7.7 % of the 8.5 % total; our measurement alone is 3.3 %. Turned around, the slope says the loop has 13.8 ± 0.3 turns.

**11 Limitation: drift (0:50).** We photographed the speeds in increasing order, with rest photos only at the start. Within each setting the phase also drifts, by up to about 8 mrad/s. A drift that is steady in time would look exactly like a slope in ω. Fitting both at once cannot separate them, because their correlation is −0.99. Next time we would take rest photos between speeds and reverse the rotation: the Sagnac phase changes sign, a drift does not.

**12 Conclusions (0:30).** The Sagnac phase is linear in ω, with a slope of 1.84 ± 0.06 that matches theory. The largest shift is half a fringe, produced by a 1 fs delay. The same slope gives c to 8.5 %, and that precision is limited by the turn count.

## Likely questions and where to go

- *Why a linear calibration when the residuals curve?* The manual asks for it. **B1** shows that a quadratic would move a and c by 1.7 %, which is inside the 3.3 % measurement uncertainty.
- *Which uncertainties are included?* **B2** lists every source (manual task 5).
- *Why leave out the random-pixel method?* **B3**.
- *What are the actual numbers?* **B4**.
- *Does the position of the rotation axis matter?* No, only the enclosed area does (Schwartz, Eur. J. Phys. 38, 015301, 2017). The loop was centred anyway.
- *Does the fibre's refractive index matter?* No. The Sagnac phase does not depend on n (same reference).
- *What about the sign of the shift?* The rings shrink, but the rotation direction was not recorded, so the sign is not compared with theory.
- *Is SMF-28e+ single-mode at 640 nm?* No. V ≈ 5.6 > 2.405, so it carries a few modes. Mode coupling that changes with stress and temperature could explain the slow phase wander. This has not been tested.

## Rebuilding

The scripts are in `source/`. From this folder, run:

```
python source/make_all.py
```

You need numpy, scipy, pandas, matplotlib and Pillow; the chaotic-circuit `.venv` has them. With ffmpeg on the PATH the GIFs come out smaller; without it, Pillow writes them. The scripts read `../sagnac_update/`, `../fringes_frames.csv`, `../calibration.json`, `../Rotation calibration.csv` and the photographs. The style (1920 × 1080, DejaVu Sans, a sentence as each headline, teal for data, orange for fits, grey for theory, small footers saying what the error bars are) copies the Chua presentation figures.
