# Handoff: Sagnac lab (lab 5), state on 23 October 2026

Written at the end of a session that was opened in the wrong project folder (`Documents/sokerproject`), so that a new conversation opened in the right place can continue. Open the new conversation with `C:\Users\danib\Documents\university-notes` as the working directory (the git repo root, whose memory folder holds the lab-5 report conventions) and point it at this file first.

## What this is

Analysis and report preparation for the Sagnac-effect experiment (Technion Physics Lab 5, measured 11 October 2026): the table's rotation is calibrated by timing turns, 161 photographs of the ring fringes at rest and nine speeds are reduced to a phase shift each, and the slope of Δφ against ω is compared with the manual's Eq. (5), Δφ = 8πNAω/(λc). Partner's report is written in two-column LaTeX in the style of the earlier Chua report (`../chaotic circuit/`).

## Where everything is

Folder: `Summer 26/lab 5/Sagnac effect/` (run every script from there, Python from `..\chaotic circuit\.venv\Scripts\python.exe`).

| File | Role |
|---|---|
| `README.md` | the data, the pipeline, the method in words |
| `RESULTS.md`, `SUMMARY.md` | every number with its uncertainty and meaning; the short version |
| `REPORT_TODO.md` | **the report plan**: the story, structure (four result sections), figures and tables column by column, section plans, abstract and conclusions plans, acknowledgments, references, traps. Numbers in it match `sagnac.txt` and `calibration.txt`. |
| `REPORT_figures_tables.tex` | the LaTeX figure and table blocks with the numbers filled in; captions left as placeholders with "must say" comments |
| `report_figures/` | the eight figures (PDF + PNG) from `report_figures.py`; seven are in the report, `rings.pdf` was dropped on 23 October |
| `PLAN_CHECKLIST.md` | the manual's tasks and preparation questions against what was done, with the report section for each |
| `uncertainty/REPORT.md` | every uncertainty component and its derivation |
| `docs/` | the code explained script by script |
| `calibration.py`, `fringes.py`, `methods.py`, `pixels.py`, `sagnac.py`, `report_figures.py` | the pipeline, in that order (`pixels.py` added 24 October); `tests/` (5 tests) |
| `previous_reports/`, `presentation_slides/`, `sagnac_update/`, `figures/` | earlier stages, kept as a record; the presentation code needs the old calibration file to run |

## Git state

Branch `sagnac-analysis-review`, five commits ahead of `main`, never merged, never pushed (the user's instruction: leave the commits on the branch). Last commit `0c1b170` "Simplify the Sagnac pipeline to library fits and plain names".

**Uncommitted** (the user was offered a commit twice and had not yet answered): the restructured report plan and its table notes. Files: `REPORT_TODO.md`, `REPORT_figures_tables.tex`, `PLAN_CHECKLIST.md`, and one sentence each in `RESULTS.md` and `uncertainty/REPORT.md` (the single-photo ring-fit precision reworded, because the simplified code no longer computes the circle-fit covariance). Suggested message: "Restructure the Sagnac report plan as one story; describe the tables column by column". Commit only when the user asks.

## Headline numbers (do not take numbers from the earlier stages)

- Calibration line ω = kS + ω₀: k = 0.02054 ± 0.00006 rad s⁻¹ per unit S, ω₀ = −0.0227 ± 0.0004 rad/s; χ²/dof = 20.2/7 = 2.9, p = 0.005 (the table bends: 3.0 % below the line at S = 90); ω at a photo setting ± 1.3 % (line's own 0.2–0.3 % ⊕ 1.28 % rms scatter of the timed settings). Reaction time 0.25 s per press since 24 October. Top speed 1.83 rad/s.
- Slope a = 1.8258 ± 0.0236 (fit) ± 0.0539 (method, half range of ring minima 1.8258 / cosine fit 1.8626 / template matching 1.7549) ± 0.0239 (calibration) = **1.83 ± 0.06 rad/(rad s⁻¹)** (3.5 %). Prediction a_th = 1.738 ± 0.137 → 1.7 ± 0.1 (N 7.7 %, R 1.8 %). a/a_th = 1.05 ± 0.09, z = 0.58, p = 0.56. N read back: 13.7 ± 0.5.
- Largest shift 3.28 ± 0.06 rad = 0.52 fringe, Δt = 1.12 ± 0.02 fs. c = (2.9 ± 0.2) × 10⁸ m/s (8.6 %), c/c₀ = 0.95 ± 0.08: the same comparison as a/a_th.
- Drift check: a = 2.1 ± 0.2 with D = −0.8 ± 0.5 mrad/s, correlation −0.99; a steady drift is not separable (speeds taken in increasing order).
- Fixed pixels (`pixels.py`, 24 October; the instructor asked for this readout, so it is in 2.3 and Fig. 7): main draw of ten pixels, slope 1.7734 ± 0.0367 (fit), s = 0.069 rad, 0 to −0.12 rad below the minima per setting; over 1001 random draws of ten, mean 1.689, SD 0.193 (11 %), 6 failed; 20/40/80 pixels: SD 0.090/0.072/0.055 (5.3/4.2/3.2 %). Shown with the SD between draws as its uncertainty, not in the method range (2.9 % stays).
- Calibration bend (`bend` block of `sagnac.txt`, added 24 October): a curve through the timed points instead of the line raises a by 2.2 % (weighted quadratic) to 2.5 % (straight segments), 0.04 rad/(rad s⁻¹); an ω² term in Δφ(ω) cannot detect it (−0.04 expected against +0.02 ± 0.05 fitted), so the ω² check is not in 2.3. Every document said "lower" until 24 October: a sign error, corrected everywhere.

## How the user wants this judged and written (feedback received, in order of weight)

1. Fit the experiment manual and the entry sheet first, then be scientifically correct, then no more complexity than needed. Over-engineered uncertainty terms ("ε ω", χ²-set scatter terms, Monte Carlo redraws, covariance scaling) were removed on request; standard methods only (SEM, reaction time of stopwatch presses, scipy `curve_fit` and `linregress`, prediction uncertainty of a line, half-range of methods).
2. Figures: no residual panels (the 19 October presentation was "despised" for them); points, fitted line and ±1σ band only; no fit results inside figures.
3. Report-facing documents carry physics and findings, not analysis mechanics ("how the sausage is made is too far"): no routine names, pixel counts, smoothing widths. Refined 24 October: a non-obvious analysis step (how the minima are found and the circles fitted, the ring numbering) is explained in a sentence or two so that nothing is a black box; routine parameters and software details stay out.
4. The report must tell one story, not list results task by task: prediction → the two measured inputs (ω in 2.1, Δφ in 2.2) → the test with its full uncertainty assembled in one section (2.3) → what limits it (2.4) → conclusions. Each section ends with the quantity the next uses; each number in one form (figure or table, not both); material that would only pad the report is listed as "left out on purpose" in `REPORT_TODO.md`.
5. Code should read like an advanced undergraduate's: library fits, pandas tables, meaningful names, nothing the result does not use (AIA, Fourier and robustness scans were deleted on 23 October).
6. Course conventions (from the Chua report): uncertainty to one significant figure and the value to the same digit, comparisons from unrounded values; every number with its kind (statistical, systematic, given); χ² always with dof and p; siunitx `\num{x \pm u}`; American spelling; the prose is the user's own (suggested text pasted verbatim had to be rewritten last time).
7. A number is quoted with its ± once, complete (24 October; the presentation was criticized for quoting the slope's uncertainty and then enlarging it once the other readouts were shown): every setting mean carries 0.04 rad from Table 2 on, and 2.3 quotes the slope once with its complete uncertainty and then explains its three contributions in turn (fit, readout, calibration), the other readouts being how the readout contribution is measured; no announcement of "three parts" (the user found one bombastic, 24 October evening), and the prediction's value is computed in 1.2, not in 2.3 (it is theory). The intercept (−0.05 ± 0.03 rad, 1.7σ) is explained, not called "consistent with zero", and the regression's uncertainty is its standard error with no χ² argument: "reduced χ² = 1 by construction" is the abandoned χ²-set scatter term in other words (rejected 24 October evening).

## Where the work stopped

The user has started writing the report's LaTeX and asked for Table 1 (the timed settings, section 2.1) to be filled in; the filled block is in `REPORT_figures_tables.tex` (Table 1) and was given in chat with two notes: 12.59 is a photo setting, not a timed one, so it does not belong in that table; and the column head N clashes with the fiber turn count N, so use n.

Likely next requests: filling Table 2 (per setting, in the tex; the budget table was dropped on 23 October), reviewing report drafts against `REPORT_TODO.md` (put each review at the top of that file, newest first, as in the Chua plan), captions, the acknowledgments' AI disclosure, and the pending commit.

Later on 23 October (evening session): the draft's 2.1 was reviewed in chat (not yet recorded in `REPORT_TODO.md`); Fig. 4 (`rings.pdf`) was dropped from the plan and the tex (block kept commented out; the plan's figure numbers stay as labels, so its Fig. 5–8 are the report's Figs. 4–7); attribution bullets were added to 1.1.2 and 2.1–2.4 of `REPORT_TODO.md`, tracing every method to the manual, Schwartz, the QE3 slides and the preparation sheet and marking the few that are ours; the two checks behind them (ring gaps in r² over 126 photos; the manual's Eq. 7 in r against r²) are recorded in `RESULTS.md`, Task 2.

On 24 October the reaction time was set to 0.25 s per press, the user's own value from an online reaction-time test (stated in draft (10)), (`REACTION_TIME_S` in `calibration.py`) at the user's request, and calibration.py, sagnac.py and report_figures.py were rerun (tests pass; figures regenerated). The headline numbers above are the new ones. In report form only Table 1's uncertainties, k, ω₀, χ², p and template matching's slope (1.75 ± 0.02) changed; a = 1.83 ± 0.06, a_th, c and Table 2 did not. Draft (5) onward carries them.

Also on 24 October, at the user's request, the uncertainty of a setting mean is presented from Table 2 on as s = 0.04 rad for every setting (the scatter of the settings about the line of 2.3, a forward reference), not the photos' standard error: `sagnac.py` writes it as `u_dphi` in `sagnac_speeds.csv`, Fig. 6's error bars show it, the tex Table 2 is regenerated from the exact means, and the plan, `RESULTS.md` and `uncertainty/REPORT.md` say so. The slope and everything after it are unchanged. The user's entry report (`Documents/Lab_5___Sagnac_Entry (6).pdf`) is reference 8 of the plan and is cited in its attribution bullets.

Drafts (3)–(10) were reviewed on 24 October; the reviews and the final checklist, with the items remaining for draft (9) at its top, are at the top of `REPORT_TODO.md`. In draft (10) the Section 2 intro, 2.1 and 2.2 are written (Tables 1–2, Figs. 1–3) with four small items open; 2.3 is next, its paragraph-by-paragraph content in the plan's checklist; 2.4, Section 1, the abstract and the conclusions are not written. The Yakobovich–Peck report (score 90) was compared on 24 October (its own section in `REPORT_TODO.md`; three small adoptions, in the checklist).

Later on 24 October: the instructor had specifically asked for the fixed-pixel (random points) readout of prep Q5, so it was ported from the presentation stage into the pipeline as `pixels.py` (same fit, current calibration; 1000 draws of ten pixels plus 60 each of 20, 40 and 80; `docs/pixels.md`; `tests/test_pixels.py`). `sagnac.py` fits its main draw like the other readouts (`fixed_pixels` in `sagnac.json`, outside the method range) and Fig. 7 shows it in both panels with the spread between draws as its uncertainty; `random_pixels.pdf` is no longer generated. Draft (11) begins 2.3; the final 2.3 content in the plan has four readouts, with the fixed pixels in their own paragraph after Fig. 7b (where their bar is the largest) and the sum in a short closing paragraph.

Draft (12) (2.3 written through the sum) was read three times late on 24 October, the third time in a fresh session, which added: the fixed pixels read about 6 % low on average (only the main draw lies inside the method range), Fig. 1 and Table 1 uncited, no citations anywhere in Section 2 yet, about twenty sentences close to the plan's suggested wording, and an AI note trimmed from the Chua one that lost its suggested-wording clause. All items are at the top of `REPORT_TODO.md`.

## Gotchas

- `python -u` for long runs; `fringes.py` takes a few minutes, everything else seconds. After a code change, rerunning the pipeline must reproduce `sagnac.json` byte for byte (it did on 23 October; on 24 October `u_dphi` and the `bend` block were added, everything else unchanged).
- `Rotation calibration.csv` starts with a BOM (`encoding='utf-8-sig'`); it was completed on 23 October (10 turns at the four fastest settings).
- The user's ω values for the photos come only from the current line; `previous_reports/`, `figures/` and `presentation_slides/source/numbers.json` quote older values (1.84 ± 0.07, top speed 1.78 or 1.81).
- Git prints CRLF warnings on commit; harmless.
