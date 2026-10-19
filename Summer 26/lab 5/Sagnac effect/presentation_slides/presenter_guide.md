# Sagnac effect: presenter's guide

This guide is for whoever presents the deck in this folder. It follows the slides in order and explains everything on them: the physics, where every number comes from, what to say, and the questions you are likely to get. It goes deeper than `presenter_notes.md`, which stays as the short version to keep open while presenting.

**How to use it.**

1. Read Part 1 (background) once. It gives you the physics and the statistics that every later slide relies on, in plain words.
2. Read Part 2 (slide by slide) before rehearsing. For each slide it lists what is on the slide, a script, the explanation in depth, where each number comes from, likely questions with answers, and the transition to the next slide.
3. On the day, keep Part 4 (numbers) and Part 3 (question bank) at hand.

**Conventions used in this guide and on the slides.**

- All uncertainties are standard uncertainties (one standard deviation, "1σ").
- On the slides an uncertainty is rounded to its first significant digit and the value to the same digit: 1.7 ± 0.1, not 1.74 ± 0.14. A number like 0.02039(4) means 0.02039 ± 0.00004.
- Intermediate numbers used in calculations, such as the regression's per-point uncertainty s = 0.045 rad and the parts of the final uncertainty (0.024, 0.054, 0.032, 0.067), are quoted unrounded.
- Every comparison (σ apart, p value) is computed from the **unrounded** numbers. This guide gives them where you need them. This matters: 1.84 and 1.7 look 0.14 apart, but the unrounded values, 1.839 and 1.738, are 0.10 apart.
- Every fit reports χ²/dof and p. A plain linear regression, whose points have no error bars of their own, reports R² instead (Part 1.5 explains why).
- The slides never cite the lab manual or the papers by name. This guide does, in brackets, so you know where each statement comes from:
  - **[Manual]** the lab manual, `Documents/Sagnac_manual.pdf`. Equation numbers as printed. The setup figure and the minimum tasks are on p. 12, the fiber-coupling section on p. 10, the fiber datasheet on p. 11.
  - **[Prep]** the preparation sheet, `Documents/Entry - Sagnac effect experiment.pdf`, questions Q1–Q12.
  - **[QE3]** the Gaussian-beam lecture slides, `Documents/QE3.pdf`.
  - **[Schwartz]** E. Schwartz, "Sagnac effect in an off-center rotating ring frame of reference", Eur. J. Phys. 38, 015301 (2017).
  - **[Pascoli]** G. Pascoli, "The Sagnac effect and its interpretation by Paul Langevin", C. R. Physique 18, 563 (2017).
  - **[Tartaglia]** A. Tartaglia and M. L. Ruggiero, "The Sagnac effect and pure geometry", Am. J. Phys. 83, 427 (2015).
  - **[Analysis]** our analysis files: `fringes.py`, `calibration.py`, `sagnac_update/`, and the presentation code in `presentation_slides/source/`. Every number in this guide is in `presentation_slides/source/numbers.json`, made by `collect_numbers.py`.

---

## The talk at a glance

**The story in five sentences.** Light sent both ways around a rotating loop takes different times, because the beam splitter moves while the light travels (slide 2). We built this with a 13-turn fiber loop on a rotating table, calibrated the table's speed with a stopwatch, and photographed the fringe pattern at ten speeds (slides 3–5). The fringes are rings, evenly spaced in r², and the phase at their center follows from a straight-line fit to the squared ring radii (slides 6–7). That phase grows linearly with ω with a slope of 1.84 ± 0.07 rad/(rad/s), 0.7σ from the predicted 1.7 ± 0.1, and the result does not depend much on how the phase is read (slides 8–10). Read through the Sagnac formula, the slope gives the speed of light, (2.8 ± 0.2) × 10⁸ m/s, limited by the given turn count (slides 11–12).

**The three numbers to remember.**

| | Value | Unrounded |
|---|---|---|
| Measured slope a | 1.84 ± 0.07 rad/(rad/s) | 1.839 ± 0.067 |
| Predicted slope a_th | 1.7 ± 0.1 rad/(rad/s) | 1.738 ± 0.137 |
| Speed of light c | (2.8 ± 0.2) × 10⁸ m/s | 2.834 ± 0.246 |

Both comparisons are 0.7σ (p ≈ 0.5): agreement.

**Timing** (about 9 minutes of talking with slide 10, 8½ without, plus transitions and questions).

| Slide | Title | Time |
|---|---|---|
| 01 | The Sagnac effect (title) | 0:20 |
| 02 | Rotation delays the beam that travels with it | 1:00 (let the animation run once) |
| 03 | Both beams travel the same 13-turn fiber loop | 0:45 |
| 04 | A stopwatch calibration converts motor setting to ω | 0:50 |
| 05 | The rings move inward as the speed is stepped up | 0:40 (let the animation play once) |
| 06 | The rings are evenly spaced in r², not in r | 1:00 |
| 07 | Squared ring radii give the phase at the center | 0:50 |
| 08 | The phase grows linearly with ω, as Sagnac predicts | 1:00 |
| 09 | Three full-image methods differ by at most 6 % | 0:50 |
| 10 | The slope from ten random pixels depends on the draw | 0:45, optional |
| 11 | The slope gives the speed of light, 0.7σ from c₀ | 0:45 |
| 12 | The Sagnac phase is linear in ω and matches theory | 0:30 |
| 13 | Thank you | questions |
| B1–B4 | Backup slides | only if asked |

**Before you present.**

- Insert the slides and GIFs as described in `README.md`. Mark B1–B4 as skipped in Google Slides (right click, Skip slide), so the talk ends on slide 13 but you can still jump to them.
- The GIFs on slides 2 and 5 start when the slide appears. Let each run once before you talk over it.
- Decide beforehand whether you show slide 10. Slide 9 already mentions the pixel method, so the talk works without it.

**How to handle questions in general.**

- Most questions fall into four groups: the physics (Part 1.1), how the phase is read (slides 6, 7, 9), the uncertainties (slides 8, 9, 11, B2) and the drift (slide 12, B4). Part 3 has short answers, sorted the same way.
- When a question is about agreement, quote the unrounded numbers and the σ.
- When you do not know, say what the data can and cannot tell. The honest answer to "is there a drift?" is "we cannot tell with this order of measurements, and here is how we would test it".

---

# Part 1. Background you need

## 1.1 The Sagnac effect

**The setup in one sentence.** Light from one source is split in two, the two beams go around the same closed loop in opposite directions and are brought together again. If the loop rotates, the two beams arrive at slightly different times, and the interference pattern shifts.

**Derivation (lab frame).** Take a circular loop of radius $R$ turning at angular velocity $\omega$, so the beam splitter on it moves at $v = \omega R$. Both beams move at $c$ in the lab frame. [Pascoli Eq. 1, Manual "classical derivation", Tartaglia Eq. 1]

- The beam going **with** the rotation has to chase the splitter. In its travel time $t_+$ the splitter moves on by $v t_+$, so the beam covers $l + v t_+ = c\,t_+$, where $l = 2\pi R$ is the loop length:
  $$t_+ = \frac{l}{c - v}$$
- The beam going **against** the rotation meets the splitter early:
  $$t_- = \frac{l}{c + v}$$
- The difference is
  $$\Delta t = t_+ - t_- = \frac{2 l v}{c^2 - v^2} \approx \frac{2 l v}{c^2} = \frac{2 (2\pi R)(\omega R)}{c^2} = \frac{4 A \omega}{c^2}, \qquad A = \pi R^2$$
  The approximation drops terms of order $(v/c)^2$.
- With $N$ turns the light goes around $N$ times, each enclosing the area $A$, and the delay is $N$ times larger [Manual Eq. 1–2 for one turn, Eq. 5 with N, Schwartz Eq. 2–3]:
  $$\Delta t = \frac{4 N A \omega}{c^2}$$
- The camera does not see a time but a phase. A delay $\Delta t$ is a phase difference of $2\pi f \Delta t$, with the light's frequency $f = c/\lambda$ [Manual Eq. 3–5, Schwartz Eq. 1]:
  $$\Delta\phi = \frac{2\pi c}{\lambda}\,\Delta t = \frac{8\pi N A}{\lambda c}\,\omega$$

This is the formula we test. It says the phase is **proportional to ω**, with the slope
$$a_{\mathrm{th}} = \frac{8\pi N A}{\lambda c} = \frac{8\pi^2 N R^2}{\lambda c}$$

**Our numbers.** $N = 13 \pm 1$, $R = 0.570 \pm 0.005$ m, $\lambda = 640$ nm, so $A = 1.02\ \mathrm{m^2}$ and

$$a_{\mathrm{th}} = \frac{8\pi^2 \cdot 13 \cdot (0.570\ \mathrm{m})^2}{640\ \mathrm{nm} \cdot 2.998\times10^8\ \mathrm{m/s}} = 1.738\ \mathrm{s} \quad\to\quad 1.7 \pm 0.1\ \mathrm{rad/(rad/s)}$$

- The unit rad/(rad/s) is "radians of phase per rad/s of rotation". It equals seconds, because the slope is really a time (it is $2\pi c/\lambda$ times $4NA/c^2$).
- **Its uncertainty.** Only $N$ and $R$ carry uncertainties. $a_{\mathrm{th}} \propto N R^2$, so the relative uncertainties add in quadrature, with $R$ counted twice because it is squared: $\sqrt{(1/13)^2 + (2 \cdot 0.5/57.0)^2} = \sqrt{7.7^2 + 1.8^2}\ \% = 7.9\ \%$, and 7.9 % of 1.738 is 0.137. λ was given without an uncertainty and $c$ is exact by definition.
- At our top speed, ω = 1.81 rad/s: $\Delta t = 4 \cdot 13 \cdot 1.02 \cdot 1.81 / c^2 = 1.07$ fs. One optical period at 640 nm is $\lambda/c = 2.13$ fs. So the delay is half a period, and the predicted phase is $3.15$ rad, close to π: half a fringe. We measured 3.29 ± 0.06 rad.
- The rim moves at only $v = \omega R \approx 1$ m/s, $v/c \approx 3 \times 10^{-9}$. The effect is still measurable because interference compares the two beams to a fraction of one optical period.

**Things the formula does and does not depend on.**

- **Only the enclosed area counts**, not the shape of the loop and not where the rotation axis is. Schwartz shows the axis position drops out: split the velocity of each point into a rotation about the loop's center plus a constant vector. The constant vector's closed-loop integral vanishes [Schwartz Eq. 4, Prep Q8].
- **The refractive index of the fiber drops out** [Schwartz Sec. 2, Prep Q10]. In the fiber the light moves at $c/n$ relative to the glass, and the glass moves with the table. With plain Galilean addition of velocities both beams would take the same time and there would be no effect at all. Relativistic addition of velocities gives exactly the vacuum result, without $n$. Historically, Fresnel's partial "drag" of light by a moving medium, measured by Fizeau, gave the same answer as an extra assumption (our reasoning. Pascoli Sec. 1.3 and 2 mention Fresnel's drag and Fizeau's measurement). Schwartz: in a medium the calculation "involves relativistic calculations", and the phase difference is independent of the material.
- **The sign follows the rotation.** Reversing the rotation reverses Δφ. We did not record the direction, so we compare only magnitudes.

**Relativity.** The slides do not mention relativity, and they do not need to. The derivation assumes only that light moves at $c$ in both directions in the lab frame, whatever the motion of the splitter. That is Einstein's second postulate. The effect is first order in $v/c$, and at that order the classical and relativistic results agree. Relativity changes Δt only at order $(v/c)^2 \approx 10^{-17}$ for us [Pascoli Sec. 3: Langevin, "an experiment of the first order in v/c", unlike Michelson–Morley. Tartaglia: Eq. 2 holds "to lowest order in v/c"]. The manual calls it "a fully relativistic effect". History: Sagnac (1913) took it as evidence for an ether. Von Laue (1920) and Langevin (1921) derived it from relativity [Pascoli Sec. 3].

**Why it matters.** Optical gyroscopes (ring-laser and fiber gyroscopes, used since the 1970s on planes, ships and submarines) measure rotation rates with exactly this effect [Manual introduction, Tartaglia Sec. I, Schwartz Sec. 1].

## 1.2 The interferometer

[Manual Fig. 1, p. 12, drawn on slide 3]

- **Laser** (640 nm) → **polarizer P₁** → **beam splitter**. The splitter sends one beam into the fiber through lens **L₁** and the other through lens **L₂**.
- Each beam goes through the whole fiber (13 turns, about 47 m) and comes out at the other end, so the two beams travel the same path in opposite directions and meet again at the splitter.
- The recombined light goes through **polarizer P₂** and lens **L₃** to the **camera**.
- **Everything turns with the table**: laser, splitter, fiber, camera. The camera therefore sees only the phase difference between the two beams, not any motion of the pattern relative to itself.
- **Why polarizers.** Light reflected from lens surfaces, the fiber holders and the splitter reaches the camera without going around the loop, and can be stronger than the signal. Two polarizers can block most of it (Malus's law), while light that went through the fiber, whose polarization the fiber changes, still gets through [Prep Q6, Manual prep question 3].
- **Why the pattern is stable.** Both beams travel the same fiber, so slow changes of the fiber length (temperature) are common to both and cancel. Rotation is different: it is non-reciprocal, it delays one direction and advances the other. That is the signal.
- **Fiber**: Corning SMF-28e+. The manual's text says SMF-28e, its printed datasheet (p. 11) and the preparation sheet say SMF-28e+. Core radius 4.1 µm, numerical aperture 0.14.

## 1.3 Interference, rings, and why r²

**Two beams, one pattern.** Where two coherent beams overlap, the intensity is
$$I = I_1 + I_2 + 2\sqrt{I_1 I_2}\,\cos\phi$$
with $\phi$ the phase difference between them [Manual Eq. 6]. Dark where $\phi = \pi + 2\pi m$, bright where $\phi = 2\pi m$. One **fringe** is a change of φ by 2π. "0.52 fringe" means $0.52 \times 2\pi = 3.3$ rad.

**Plane waves give straight fringes.** Two plane waves meeting at a small angle give straight, evenly spaced fringes with period $\Lambda = \lambda/(2\sin\theta)$, with θ each beam's angle to the axis [Manual Eq. 6 and Fig. 2, Prep Q1–Q2]. That is the textbook picture of a Sagnac interferometer.

**Our beams are Gaussian beams.** They leave the fiber as Gaussian beams, whose wavefronts are spherical with radius $R(z)$. Such a beam carries the transverse phase $k\rho^2/(2R(z))$ at distance ρ from its axis, plus an on-axis phase $kz - \zeta(z)$ ($\zeta$ is the Gouy phase) [QE3, "The phase of the Gaussian Beam", PDF p. 10, Prep Q3].

**Two Gaussian beams with different wavefront radii give rings.** If the two beams reach the camera on the same axis but with wavefront radii $R_1 \neq R_2$, their phase difference at radius $r$ is
$$\phi(r) = \frac{k r^2}{2}\left(\frac{1}{R_1} - \frac{1}{R_2}\right) + \phi_c = 2\pi\frac{r^2}{P} + \phi_c, \qquad P = \frac{2\lambda}{|1/R_1 - 1/R_2|}$$
using $k = 2\pi/\lambda$ [Analysis, `fringes.py` header].

- $\phi_c$ is the phase difference at the center of the pattern.
- $P$ is the ring spacing **in $r^2$**: going out by $P$ in $r^2$ adds one fringe. It is the same for every ring and is set only by the two wavefront radii.
- So the dark rings sit at $r_m^2 = P(m + \eta)$: equally spaced in $r^2$, and closer and closer together in $r$ (like Newton's rings). In our fastest photo the rings are 241 px and then 132 px apart in $r$, but 1.08 × 10⁵ px² both times in $r^2$.
- Schwartz notes the same: his pattern is not linear like Young's fringes, and fringe widths differ between the center and the edge [Schwartz Sec. 3].

**What moves the pattern.**

| Change | Effect on the rings | Effect on φ_c |
|---|---|---|
| Rotation | the same phase is added everywhere, so the rings slide in or out | changes, this is the Sagnac signal |
| Relative tilt of the two beams | adds a phase linear across the camera, which mainly moves the ring center | small, only through a center shift |
| A lens-to-fiber distance changes | the wavefront radii change, so P changes | small, through the Gouy phase |

That is why the analysis tracks the center phase $\phi_c$ and fits the center and $P$ in every photo [Analysis, `fringes.py` header, Prep Q7 asks about the lens-to-fiber distance].

**Which way the rings move** depends on the sign of $1/R_1 - 1/R_2$ and on the rotation direction. In our run they move inward as ω grows.

## 1.4 Fiber modes and coupling

- A step-index fiber's **V number** is $V = (2\pi/\lambda)\,a\,\mathrm{NA}$, with $a$ the core radius [Manual Eq. 9]. For $V < 2.405$ the fiber guides a single mode, above it several.
- With $a = 4.1$ µm and NA = 0.14 the manual itself gets $V = 5.6$ at 640 nm [Manual, "Optimum fiber coupling", p. 10]. Its V-number section works the same numbers out for 532 nm (6.78) and 1550 nm (2.32). The datasheet gives the cutoff wavelength, 1260 nm: below it the fiber is not single-mode. So at 640 nm our fiber guides several modes, and Prep Q12, which asks whether SMF-28e+ is single-mode at 640 nm, is answered no.
- How much light goes into the **fundamental mode** depends on how well the focused beam matches it. The manual computes up to 93 % with a matched lens (f = 15.5 mm, waist 3.16 µm) and at most 47 % with the setup's 8 mm lens [Manual p. 10]. Coupling into one mode is possible with this fiber, but in practice it was very hard to align.
- **Why it matters.** In one mode, the two directions see exactly the same fiber, so stress and temperature cancel (reciprocity). With several modes, power moves between modes along the fiber, not necessarily the same way in both directions, and the camera records the whole multimode field. Then a slowly varying phase difference can appear that has nothing to do with rotation. Fiber gyroscopes put a single-mode filter at the common port for this reason (general knowledge, not in our documents). This is our plausible, untested explanation for the phase wander (slides 8, 12, B4).

## 1.5 Statistics toolkit

Every statistics term on the slides, in plain words.

**Mean, SD, SEM.** For $n$ measurements $x_i$: the mean $\bar x$, the standard deviation $\mathrm{SD} = \sqrt{\sum (x_i - \bar x)^2/(n-1)}$ (how much a single measurement varies), and the **standard error of the mean** $\mathrm{SEM} = \mathrm{SD}/\sqrt n$ (how well the mean is known, **if** the measurements are independent). Slide B3 and the B2 table use SEM.

**Standard uncertainty.** "± u" means one standard deviation. For a normal distribution the true value lies within ±u about 68 % of the time.

**Combining uncertainties.**

- Independent sources add **in quadrature**: $u = \sqrt{u_1^2 + u_2^2 + \dots}$. Reason: variances of independent quantities add. If the sources were fully correlated, they would add linearly (the worst case).
- For a product or ratio, the **relative** uncertainties add in quadrature, and a power counts as many times as its exponent. For $c = 8\pi^2 N R^2/(\lambda a)$:
  $$\left(\frac{u_c}{c}\right)^2 = \left(\frac{u_N}{N}\right)^2 + \left(2\frac{u_R}{R}\right)^2 + \left(\frac{u_a}{a}\right)^2$$
- **Statistical** uncertainties shrink with more data (photo noise). **Systematic** ones do not (choice of method, choice of calibration model). For a systematic spread between a few alternatives we take half the range as one standard uncertainty. That is a convention. A flat distribution over the full range would give range/√12, a cautious choice the full range.

**Linear regression** (least squares, equal weights). Fit $y = a x + b$ by minimizing $\sum (y_i - a x_i - b)^2$. This is how Schwartz gets his slopes [Schwartz Eq. 7].

- The points have no error bars of their own. The regression measures how far they lie from the line and uses that as the uncertainty of every point:
  $$s = \sqrt{\frac{\sum (y_i - \text{line}_i)^2}{n - 2}}$$
  The $n - 2$ (not $n$) is like the $n - 1$ of an SD: two numbers (slope and intercept) were already taken from the same points.
- The slope's standard error is $u_a = s / \sqrt{\sum (x_i - \bar x)^2}$. It shrinks when the points lie close to the line and when the $x$ values are spread out.
- **R²** $= 1 - \sum(y_i - \text{line}_i)^2 / \sum (y_i - \bar y)^2$: the fraction of the variation of $y$ that the line accounts for. 1 means every point on the line.

**χ² test** (a fit to points **with** their own error bars $u_i$).

- $\chi^2 = \sum (y_i - \text{fit}_i)^2 / u_i^2$. If the model is right and the error bars are right, each term is about 1 on average, so $\chi^2 \approx$ dof, where dof = number of points − number of fitted parameters. So χ²/dof ≈ 1 is good.
- **p** is the probability of getting a χ² at least this large if the model and the error bars were right. Tiny p (10⁻⁸) means the model, or the error bars, are not enough. Very large p (near 1) can mean the error bars are too big.
- **χ² or R²?** Use χ² when the points have trustworthy error bars of their own (slide 4: the stopwatch). In a plain regression the noise is taken from the points themselves, so χ²/dof = 1 by construction and says nothing. Then R² is the fit-quality number (slides 8, 9). Note that R² can be 0.999 for a line that χ² rejects: R² does not know how small the error bars are.

**Comparing two values: σ and p.**

- $z = (x_1 - x_2)/\sqrt{u_1^2 + u_2^2}$ is the difference in units of the combined uncertainty ("0.7σ apart").
- The **two-sided p** is the probability of a difference at least this large, in either direction, if the two values were really equal and the uncertainties right: $p = 2(1 - \Phi(|z|))$, with Φ the normal distribution. 0.7σ gives p ≈ 0.5, 1σ gives 0.32, 2σ gives 0.05, 3σ gives 0.003.
- p is **not** the probability that the theory is right. A large p only says the data give no reason to doubt the agreement.
- Rule of thumb: below about 2σ agree, above about 3σ disagree.

**Testing one fitted parameter (t test).** $t$ = estimate / its standard error, compared with a t distribution with dof = points − parameters. Used for the drift term D on B4 (p = 0.12), the intercept and the added ω² term on slide 8.

**Residuals.** The differences between the data and a fitted curve ("Δφ − line" on slide 8, "measured − model" on B1).

**Correlation.**

- Between two data sets: $r$ from −1 to 1, how well one follows the other linearly (used in the slide 8 checks).
- Between two fitted parameters: ρ = −0.99 (B4) means the fit cannot tell the two parameters apart: raising one and lowering the other fits almost equally well.

**Alternating least squares** (slide 10). A fit with two sets of unknowns, here the pixel parameters and the photo phases. It fixes one set and fits the other, then swaps, and repeats until nothing changes.

**Bootstrap.** To get the uncertainty of a summary number such as the SD of 1000 slopes, resample the 1000 values with replacement many times and see how much the summary varies (slide 10's (5) on 0.196).

**Robust spread.** Half the width of the middle 68 % of a set of values. For a normal distribution it equals the SD, but a few outliers hardly move it (slide 10: 0.187 against SD 0.196).

---

# Part 2. Slide by slide

Each slide has the same headings: **Purpose**, **On the slide**, **Say** (a script you can adapt), **In depth**, **Numbers and where they come from**, **Questions you may get**, **Transition**.

## Slide 01. The Sagnac effect (title)

**Purpose.** Name the experiment and show, before any explanation, that the photos change with rotation.

**On the slide.**
- Title "The Sagnac effect", subtitle "The Sagnac phase in a rotating fiber loop, and the speed of light".
- Names and instructor.
- Five real photographs from our run, at rest and at 0.40, 0.88, 1.28 and 1.81 rad/s, each centered on its fitted ring center.

**Say.** "We measured the Sagnac effect: light sent both ways around a fiber loop on a rotating table. These are five of our photographs. From left to right the table turns faster, and the dark ring closes in on the center. From the size of that shift we test the Sagnac formula and get the speed of light."

**In depth.** The five photos are the representative photos of the settings S = 0, 20.74, 44.07, 64.00 and 90. The center goes from bright (rest) to nearly dark (top speed): a phase change of about half a fringe.

**Transition.** "Why should rotation change anything? Here is the idea."

## Slide 02. Rotation delays the beam that travels with it

**Purpose.** The physics in one picture and two formulas. Let the animation run once before talking over it.

**On the slide.**
- Left, animation (loops every 7.3 s, not to scale): a turning loop with a splitter on it, two pulses leaving the splitter in opposite directions. The blue pulse goes with the rotation, the orange one against it. The bar at the bottom marks when each returns to the splitter, and the gap Δt between them.
- Right, top: "Delay after N turns, each enclosing an area A:" $\Delta t = 4NA\omega/c^2$.
- Right, middle: "Why the camera sees it as a phase shift:" two waves of the same light, orange (against rotation) and blue (with rotation), the blue one shifted later by Δt. One period is marked, $T = \lambda/c$. No text with it: the explanation is spoken (see Say).
- Right, bottom: $\Delta\phi = (2\pi/T)\,\Delta t = 8\pi N A\omega/(\lambda c)$.
- Footer: "Animation not to scale. The rim is drawn moving at 0.25c and the loop as a single turn."
- Subtitle: optical gyroscopes, delay proportional to the rate of rotation.

**Say.** "A beam splitter on the table sends light both ways around a closed loop. While the light travels, the splitter moves with the table. The beam going with the rotation, blue, has to chase the splitter, so it needs a little more than one turn. The beam going against the rotation, orange, meets the splitter early. The difference in return times is Δt, proportional to the enclosed area, the number of turns and the rate of rotation. Why does the camera see a delay? Both beams carry the same laser light, a wave with period T. Delaying one beam by Δt shifts its wave by the fraction Δt/T of a period: a phase shift of 2π Δt/T. Where the two beams recombine they form bright and dark bands, and a phase shift moves those bands. So the phase shift is Δt times 2π/T, proportional to ω."

**In depth.** The full derivation is in Part 1.1. Key points if asked:
- In the lab frame both beams move at c. The co-rotating one needs $t_+ = l/(c - v)$, the other $t_- = l/(c + v)$. $\Delta t = 2lv/(c^2 - v^2) \approx 4A\omega/c^2$ per turn.
- **Why a phase.** Both beams carry light of one frequency $f = 1/T = c/\lambda$. At the recombination point, beam 2's field is $\cos(2\pi f t)$ and the delayed beam 1's is $\cos(2\pi f (t - \Delta t))$. They differ by the phase $2\pi f\,\Delta t = 2\pi\,\Delta t/T$. The camera records their interference, brightness $\propto 1 + \cos(\text{phase difference})$, so changing the phase difference moves the bright and dark bands. How that turns into rings moving inward is slide 6.
- The size of the effect is deliberately not on this slide, because λ and fringes come later. If asked: at our top speed Δt = 1.07 fs, about half of one light period (T = 2.13 fs at 640 nm), so we expect about half a fringe. Say it on slide 5, where the rings visibly move by about half a fringe.
- The formula is linear in ω. That is the main thing we test.

**Numbers and where they come from.**
- $\Delta t = 4 \cdot 13 \cdot (\pi \cdot 0.570^2\ \mathrm{m^2}) \cdot 1.81\ \mathrm{rad/s} / c^2 = 1.07$ fs.
- $\lambda/c = 640\ \mathrm{nm}/c = 2.13$ fs.
- Predicted $\Delta\phi = 2\pi \cdot 1.07/2.13 = 3.15$ rad (shown as 3.2), close to π. We measured 3.29 ± 0.06 rad (slide 12).
- In the animation the rim moves at 0.25c, so the return times are $2\pi/(1 - 0.25)$ and $2\pi/(1 + 0.25)$ in units of R/c. Our real rim moves at about 1 m/s, $3 \times 10^{-9}\,c$.

**Questions you may get.**
- *Is this relativity?* The derivation needs only that light moves at c both ways in the lab frame, Einstein's second postulate. It is a first-order effect in v/c, where the classical and relativistic results agree. Corrections are of order $(v/c)^2 \approx 10^{-17}$ [Pascoli Sec. 3, Tartaglia]. Part 1.1 has more.
- *Does the fiber's refractive index change the formula?* No. Relativistic velocity addition makes n drop out (Part 1.1) [Schwartz Sec. 2, Prep Q10].
- *Does it matter where the axis is?* No, only the enclosed area [Schwartz Eq. 4, Prep Q8].
- *Why A of one turn times N, not the total area of the fiber?* Each turn encloses A, and the light goes around N times, so the delays add.

**Transition.** "This is how we realized it."

## Slide 03. Both beams travel the same 13-turn fiber loop

**Purpose.** The apparatus, the given parameters, and the prediction we will test.

**On the slide.**
- Left: top view of the setup, after the manual's setup figure. Laser → P₁ → beam splitter. L₁ and L₂ couple the light into the two fiber ends. The two beams (blue beam 1, orange beam 2) go around the 13-turn loop in opposite directions. P₂ and L₃ lead to the camera. A dark arrow marks ω.
- Right, "SYSTEM PARAMETERS (given)": fiber turns N = 13 ± 1, loop radius R = 57.0 ± 0.5 cm, wavelength λ = 640 nm, fiber SMF-28e+.
- "WHAT WE MEASURE": angular velocity ω with a stopwatch, phase shift Δφ from camera photos.
- "PREDICTION": $\Delta\phi = (8\pi NA/\lambda c)\,\omega \equiv a_{\mathrm{th}}\,\omega$, $a_{\mathrm{th}}$ = 1.7 ± 0.1 rad/(rad/s), with A = πR², the area of one turn.
- Footer: P = polarizer, L₁ and L₂ fiber coupling lenses, L₃ camera lens. The polarizers reduce reflections. λ was given without an uncertainty.

**Say.** "The laser goes through a polarizer to the beam splitter. One beam enters the fiber at L₁, the other at L₂. They travel the same 13 turns in opposite directions, recombine at the splitter and go to the camera. Everything you see turns with the table. We measure two things: ω with a stopwatch, and Δφ from photographs of the fringes. With the given loop parameters the prediction is a straight line with slope 1.7 ± 0.1 radians per rad/s."

**In depth.**
- The drawing follows the manual's Fig. 1 (p. 12). That figure also shows a second, small loop (R = 15.5 ± 0.1 cm, 46 ± 1 turns) used for task 4. It is left out because task 4 was skipped on the instructor's instruction.
- The polarizers are there to suppress reflections that reach the camera without going around the loop [Prep Q6, Manual prep question 3]. See Part 1.2.
- **The prediction.** With A = πR² = 1.02 m², $a_{\mathrm{th}} = 8\pi^2 N R^2/(\lambda c) = 1.738 \pm 0.137$, rounded to 1.7 ± 0.1 [Manual Eq. 5].
- **The ±0.137.** From N (1/13 = 7.7 %) and R (squared: 2 × 0.5/57.0 = 1.8 %), $\sqrt{7.7^2 + 1.8^2} = 7.9\ \%$. The slide does not explain it. If asked, say it verbally: almost all of it is the turn count (7.7 %), the radius adds 1.8 %. The same numbers reappear in the budget on slide 11.
- **The fiber name.** The manual's text says SMF-28e (V-number section, the heading of its specification page, its preparation task). The datasheet printed under that heading is Corning's SMF-28e+, and the preparation sheet says SMF-28e+ (Prep Q11–Q12). The core radius 4.1 µm and NA 0.14 come from that datasheet. Nothing in the analysis depends on the name.

**Numbers and where they come from.** N, R and λ are given in the caption of the manual's Fig. 1 (p. 12). $a_{\mathrm{th}}$ from Manual Eq. 5 with A = πR².

**Questions you may get.**
- *Where does 13 ± 1 come from? Did you count the turns?* It is given in the manual's setup figure. We did not count them. We treat ±1 as one standard deviation. If it is meant as a bound (13 for sure within ±1), the standard uncertainty would be 1/√3 = 0.58 turns, and the uncertainty of c would drop from 8.7 % to about 6 %. Keeping ±1 is the cautious choice.
- *What is R exactly?* The given loop radius, taken as the mean radius of the windings.
- *Why polarizers?* See above and Part 1.2.
- *Is the fiber single-mode?* Not at 640 nm (V = 5.6). See slide 12 and B4.

**Transition.** "To use the formula we need ω. The motor controller does not give it, so we calibrated it."

## Slide 04. A stopwatch calibration converts motor setting to ω

**Purpose.** How ω is known, with its uncertainty, and why we use a straight line although it does not fit perfectly.

**On the slide.**
- Subtitle: mean turn time T at nine settings S of the motor controller (S sets the motor voltage), ω = 2π/T.
- Left: timed ω against the motor setting S (blue points), the linear fit (orange), and dashed vertical lines at the nine settings used for the photos.
- Right, "LINEAR FIT": $\omega = kS + \omega_0$, k = 0.02039(4) rad/s per unit S, χ²/dof = 47.2/7, p = 5 × 10⁻⁸, "points within 2.3 % of the line, error bars 0.5 % of ω".
- "ONE POINT: S = 90, 8 TURNS TIMED", "Each dot is one turn (stopwatch splits)": the eight turn times, their mean ⟨T⟩ = 3.545 s (orange line) and ω = 2π/⟨T⟩ = 1.772 rad/s.
- The error bars of the timed points (stopwatch 0.04–0.18 %, reproducibility ε 0.48 %, explained below) are not listed on the slide. They are smaller than the markers.
- Footer: the photos take ω from this line: 1.81 rad/s at S = 90, 2.3 % above the timed 1.772 rad/s. This difference is part of the final uncertainty (calibration model).

**Say.** "The motor controller gives a setting S, not an angular velocity. At nine settings we timed five to eight full turns with a stopwatch, turned the mean turn time into ω = 2π/T, and fitted a straight line, as the first task asks. The line does not fit within the tiny error bars: χ²/dof is 47 for 7, because the table turns a little slower than linear at high settings. We still use the line, and on slide 9 we include the difference a curved calibration would make in our uncertainty."

**In depth.**
- **Task 1** asks for the (linear) ratio between the motor voltage and ω [Manual minimum task 1]. S is the number set on the motor controller (no unit), which sets the motor voltage. We calibrated ω directly against S.
- **Stopwatch splits.** The turns were taken as splits: one press of "lap" at each pass, so neighboring turns share a press. The mean turn time is $T = (t_n - t_0)/n$, so only the first and last press matter. If each press has an independent error $\sigma_r$, then $\mathrm{Var}(T) = 2\sigma_r^2/n^2$, while the individual turn times have $E[\mathrm{SD}^2] = 2\sigma_r^2(n+1)/n$. Eliminating $\sigma_r$:
  $$u(T) = \frac{\mathrm{SD}}{\sqrt{n(n+1)}}$$
  [Analysis, `calibration.py` header]. Check that the turns really are splits: neighboring turn times should be anticorrelated (a late press makes one turn long and the next short). Measured −0.34, expected −0.43 ± 0.09 for splits.
- Since ω = 2π/T, $u(\omega)/\omega = u(T)/T$. At S = 90: SD = 0.054 s, u(T) = 0.054/√72 = 0.006 s, 0.18 %. Over the nine settings 0.04–0.18 %.
- **Speed reproducibility ε = 0.48 %.** Set the table to the same S again and its real speed differs by about 0.48 %. This is the spread of the settings about a smooth curve that the stopwatch cannot explain. It is set by χ²/dof = 1 for a quadratic fit, and added in quadrature to every error bar and later to every ω of a photo.
- **Why χ², not R², here.** These points have error bars of their own, so χ² tests the line against them, the stronger test. R² would be above 0.999 for this rejected line. (Part 1.5.)
- **How to defend the χ².** The line is rejected because the timing is very precise, not because it is far off. The points lie within 2.3 % of the line (−0.6, +1.3, +0.8, +0.7, +1.1, +0.7, −0.4, −1.3, −2.3 % in order of S), while each error bar is only about 0.5 %. The deviations form a bend (above the line in the middle, below at the top), which a quadratic removes (χ²/dof = 6.0/6). Say: "The stopwatch is precise enough to see that the table turns about 2 % slower than linear at the top. We use the line, as the task asks, and carry what a curved calibration would change, 1.7 % on the slope, as an uncertainty. It cannot change the agreement." Do not switch to R² here: R² = 0.9991 for this rejected line, because R² ignores the error bars.
- **Why the line anyway.** The task asks for the linear relation. The rejection is about 2 % effects: at S = 90 the line gives 1.81 rad/s, the timed value is 1.772. The photos take ω from the line. A quadratic follows the points (B1) and raises the Sagnac slope by 1.7 %, which we carry as the "calibration model" uncertainty (slide 9).
- The intercept $\omega_0 = -0.0212(6)$ rad/s is not shown. The line was only tested inside S = 5.78–90, so its value at S = 0 means nothing. At rest ω = 0 exactly.

**Numbers and where they come from.** `Rotation calibration.csv` (turn times), `calibration.json` and `calibration.txt` (fits). Stopwatch-only line: χ²/dof = 807/7. With ε: 47.2/7, p = 5 × 10⁻⁸. Quadratic with ε: 6.0/6 by construction. Cubic: 2.6/5, p = 0.75.

**Questions you may get.**
- *Why not use the quadratic, since it fits?* The task asks for the linear relation, and the quadratic is fitted to the same nine points it is then judged on (its χ²/dof = 1 is by construction). We show both (B1) and include the 1.7 % difference in the uncertainty.
- *Why are the error bars invisible?* At most 0.5 % of ω, smaller than the markers.
- *Was ω measured during the photos?* No. Each photo gets the calibrated ω of its setting. ε covers the variation between visits to the same setting.
- *How do you know ε?* From the calibration itself: it is the extra spread needed for a smooth curve to fit.

**Transition.** "With ω known, here is everything we photographed."

## Slide 05. The rings move inward as the speed is stepped up

**Purpose.** Show the raw data before any analysis. Let the animation play once.

**On the slide.**
- Subtitle: all 161 photographs in recording order, at rest, then nine speeds in increasing order.
- Animation (loops every 19.4 s, 0.1 s per photo, first photo held 1 s, last 2.5 s). Top row: ω of the current photo, motor setting S, "photo j of n" within the setting, and the camera clock.
- Left: the photo, cropped on the ring center at rest, with a dashed circle at the first dark ring at rest.
- Right: "Rate of rotation during the run", ω against minutes since the first photograph: a staircase that only goes up. A marker follows the animation.
- Footer: ω of each photo from the calibration on the previous slide. Dashed circle: the first dark ring at rest. Photos lightly smoothed for display.

**Say.** "Before any analysis, this is all of our data: 161 photographs in the order we took them. The dashed circle is where the first dark ring was at rest. As the table speeds up, the rings are pulled in toward the center, by about half a ring spacing at the top speed. That is the size slide 2 predicts: a delay of about one femtosecond, half a period of our 640 nm light. Why they move inward, rather than just getting brighter or darker, is on the next slide. On the right is the rate of rotation during the run: ten levels in twelve minutes, always increasing."

**In depth.**
- 10 photos at rest and 15–22 at each of nine speeds, 16:21:40 to 16:33:46.
- The first dark ring shrinks from 258 px at rest to 101 px at the top speed. It is the same ring moving in, not a new one appearing, so the phase changes by less than one fringe (about half).
- The title says "as the speed is stepped up", not "because of rotation". Speed and elapsed time increase together in this run, which is the limitation of slide 12.
- The crop is fixed, so you can also see the ring center wander by a few pixels (3 px SD over the run). That is a small change of the relative tilt of the two beams. It is why the analysis fits the center in every photo.
- Display only: the photos are lightly smoothed. The analysis uses the original photos. The camera never saturates (maximum 239 of 255).

**Questions you may get.**
- *Why only half a fringe?* That is what the formula predicts for our loop and speeds: Δt = 1.07 fs at 1.81 rad/s, against a light period of T = λ/c = 2.13 fs, so a phase of 2π × 1.07/2.13 = 3.15 rad.
- *Why increasing order only?* That is how the run was done. It is our main limitation (slide 12).

**Transition.** "To measure how far the rings moved, we first need to know where they are."

## Slide 06. The rings are evenly spaced in r², not in r

**Purpose.** How the rings are found, how good the fit is, and the key observation that they are evenly spaced in r².

**On the slide.**
- Subtitle: fastest-rotation photo, dark-ring minima found along 170 rays, concentric circles fitted to them.
- Left: the photo with the rays, the colored dots (the lowest point of each dark ring along each ray, one color per ring), the white dashed fitted circles labeled m = 0, 1, 2, the center (+), and "circle fit, rms 6 px".
- Under it: "Two beams with different wavefront radii:" $\phi(r) = 2\pi r^2/P + \phi_c$. "Rotation shifts φ_c, so the rings move. P is the ring spacing in r², the same for every ring."
- Right, top: intensity against r. The dark rings are 241 px and then 132 px apart: "the spacing shrinks outward".
- Right, bottom: the same intensity against r². Both spacings are 1.08 × 10⁵ px²: "one constant spacing P".
- Footer: rays every 1.5° from −75° to 180° (no fringes at lower left). Dots: darkest point of each ring on each ray. All photos: rms 3–8 px, radii typically ±1 px. $P = 2\lambda/|1/R_1 - 1/R_2|$ for wavefront radii R₁, R₂. Median over all photographs: P = 1.07 × 10⁵ px².

**Say.** "To measure the rings in every photo, we send 170 rays out from the center and find where each dark ring is lowest along each ray: these are the colored dots. Then we fit concentric circles to all the dots at once, which gives the center and the radius of every ring. The dots lie within a few pixels of their circles. And the radii hold a surprise: the rings are not evenly spaced in r. Against r² they are, with one constant spacing P. That is what two Gaussian beams with different wavefront radii give: their phase difference grows in proportion to r², and the dark rings sit where it reaches π, 3π, 5π. And this is why rotation moves them inward: rotation adds the same phase Δφ everywhere. To stay at π, 3π, 5π, every dark ring has to move to a smaller r². The pattern keeps its shape and slides inward, and after one full fringe each ring would sit where its inner neighbor was."

**In depth.**

*The ring finder* [Analysis, `fringes.py` steps 1–4]:
1. The red channel (the only one with light) is smoothed with a Gaussian of σ = 4 px against speckle.
2. 170 rays every 1.5° from −75° to 180°. The lower-left part of the pattern has no fringes.
3. Along each ray, each dark ring's minimum is searched within ±0.3 of the local fringe period around a first guess, then refined to below a pixel with a parabola through the neighboring samples. Shallow minima (noise) are dropped.
4. Concentric circles (one center, one radius per ring) are fitted to all the minima with a robust loss, so a few stray points have little pull. The rays are re-centered on the new center and the search is repeated, three rounds in all.

*How good the fit is.*
- "rms 6 px" is the root-mean-square distance of the single dots from their circles in this photo. Over all 161 photos 2.8–8.2 px, median 4.7 px. This is speckle on single points, small next to the 130–240 px between rings.
- "radii typically ±1 px" is the uncertainty of each fitted radius. Each ring has 54–170 minima (median 136), so its radius is known much better than any single dot: 0.3–2.3 px, median 0.8 px. The center is known to about ±0.6–0.7 px.
- No χ² here: the single minima have no independent error bars, so their measured distance from the circles is their noise.

*Why rings, and why r²* (Part 1.3 has the full argument):
- Two plane waves give straight fringes. Our beams leave the fiber as Gaussian beams with spherical wavefronts, phase $k\rho^2/(2R(z))$ [QE3, PDF p. 10, Prep Q3].
- With wavefront radii R₁ ≠ R₂ at the camera, the phase difference is $\phi(r) = 2\pi r^2/P + \phi_c$ with $P = 2\lambda/|1/R_1 - 1/R_2|$.
- Going out by P in r² adds one fringe, so the dark rings are equally spaced in r². In r they crowd together outward.
- Schwartz noticed the same non-uniform spacing [Schwartz Sec. 3].
- **Why the rings move radially** (not on the slide, say it). Rotation delays the whole beam, so it adds the same Δφ at every point: φ_c → φ_c + Δφ. A dark ring needs the total phase to stay at π + 2πm, so its r² must drop by exactly the added amount: $R_m^2 \to R_m^2 - P\,\Delta\phi/2\pi$, the same step for every ring. Analogy: the rings are contour lines of a bowl, and rotation fills the bowl evenly, so every contour shrinks toward the center and none appears or disappears. Our half fringe (Δφ = 3.3 rad) moves the first dark ring from 258 px to 101 px.
- Rotation changes only φ_c. A tilt between the beams mainly moves the center. A change of a lens-to-fiber distance mainly changes P [Prep Q7 asks about this]. Neither is completely free of phase: a tilt about a pivot away from the center shifts φ_c a little, and refocusing changes the Gouy phase on the axis [QE3]. Such effects are a plausible part of the wander between settings (slide 8).
- The profiles on the right are the intensity averaged over the fringe sector in thin rings around the fitted center.

**Numbers and where they come from.** `fringes_frames.csv` (per-photo center, radii, their uncertainties, rms, P). Spacings measured on this photo's profile.

**Questions you may get.**
- *Why not fit straight fringes, like the manual's Eq. 6?* Because the photos show rings, and the ring spacing follows r², not r. The straight-fringe formula is the plane-wave limit.
- *What sets P? Can you get R₁ and R₂ from it?* The difference 1/R₁ − 1/R₂ of the two wavefronts. To get the radii in meters you would need the pixel size and the magnification, which the phase does not need (slide 7).
- *Why is the lower left missing?* Our photos have no fringes there, so the rays leave it out.
- *Why a robust fit?* Speckle sometimes produces a false minimum. A robust loss lets those few points pull less than ordinary least squares would.

**Transition.** "If the rings are equally spaced in r², their radii tell us the phase at the center directly."

## Slide 07. Squared ring radii give the phase at the center

**Purpose.** How one number, the center phase, is extracted from each photo.

**On the slide.**
- Subtitle: one photograph at rest and one at the highest speed, the same fit applied to all 161.
- Left: the two photos with their fitted dark rings (blue frame: at rest, orange frame: ω = 1.81 rad/s).
- Middle: $R_m^2/P$ against the dark-ring number m, for both photos. Two parallel straight lines. At m = 0 they are 0.52 fringe apart.
- Right: "FIT EACH PHOTO" $R_m^2 = P(m + \eta)$, slope P, intercept Pη. "PHASE AT THE CENTER" $\phi_c = \pi - 2\pi\eta$, a dark ring has phase π. "THIS PAIR" η: 0.62 → 0.10, $\Delta\phi = 2\pi \times 0.52 \approx 3.3$ rad.
- Footer: per photo, center and radii from the previous slide (rms 3–8 px over all photos), then a straight line in $R_m^2$. P is fitted for every photo, so a change of ring spacing does not move φ_c.

**Say.** "For every photo we now have the ring radii. Plotting R² against the ring number gives a straight line, as we just saw: the slope is the spacing P, and the intercept tells us how far the innermost dark ring is from the center, in units of one fringe. That fraction, η, is the phase at the center. Here η goes from 0.62 at rest to 0.10 at the top speed: a shift of 0.52 fringe, 3.3 radians."

**In depth.**
- From $\phi(r) = 2\pi r^2/P + \phi_c$, a dark ring is where $\phi = \pi + 2\pi m$ [Manual Eq. 6 with cos = −1, Prep Q1]. Solving for r²:
  $$R_m^2 = P\left(m + \eta\right), \qquad \eta = \frac12 - \frac{\phi_c}{2\pi}, \qquad \phi_c = \pi - 2\pi\eta$$
- In practice: fit $R_m^2 = Pm + b$ for each photo, then $\eta = b/P$. The plot shows $R_m^2/P$, so the slope is 1 by construction and the intercept is η.
- With two or three rings per photo this line has at most one degree of freedom, so no fit statistic would mean much. The quality check is the circle fit of slide 6, which uses all the minima.
- $\Delta\phi = \phi_c - \phi_c(\text{rest})$, with the rest value the mean of the 10 rest photos. For this pair: $2\pi(0.62 - 0.10) = 2\pi \times 0.52 = 3.3$ rad.
- This is the manual's method (shift of each minimum over the local fringe period, averaged over the minima) [Manual Eq. 7] and Schwartz's $\Delta\phi = 2\pi\,\Delta x/\delta$ [Schwartz Eq. 5], done in r², where every ring has the same period, so one line uses all the rings at once.
- **Why it is robust.** P is fitted in every photo, so a change of the ring spacing does not move φ_c. Minimum positions do not depend on the camera's brightness response.
- **No pixel-to-length calibration is needed.** η = b/P is a ratio of two numbers in px². If a pixel were a length s, both would be multiplied by s² and η would not change. A length calibration would only be needed to turn P into the physical wavefront radii.
- **Ring numbering.** η stays between 0 and 1 for the whole run (0.62 at rest, 0.10 at the top speed), so the numbering never has to shift. At rest only two rings are in the fringe region, so the line goes exactly through both points there. At the three lowest speeds some photos have two, from S = 44 on every photo has three.

**Numbers and where they come from.** η and P per photo: `fringes_frames.csv`. The two photos are the representative photos of S = 0 and S = 90.

**Questions you may get.**
- *With only two rings at rest, how do you know the line is straight?* You do not, from that photo alone. Equal spacing in r² is checked on the photos with three rings (slide 6) and is what the beam model predicts.
- *What if a ring crosses the center?* Then η would wrap from 0 to 1 and the numbering would shift by one. It does not happen in our run.
- *Why the center phase and not the phase at some radius?* Rotation changes only φ_c, while a change of P changes the phase at every other radius. The center is the one place where a change of ring spacing adds nothing.

**Transition.** "Doing this for all 161 photos gives one phase per photo. Averaged per setting, here they are against ω."

## Slide 08. The phase grows linearly with ω, as Sagnac predicts

**Purpose.** The main result: a straight line, its slope, and the comparison with the prediction.

**On the slide.**
- Upper left: the ten setting means of Δφ against ω (blue), the linear regression (orange), and the prediction ± 1σ (gray dashed line and band). The right axis gives the same phase in fringes (Δφ/2π).
- Lower left: each point minus the line, "Δφ − line (rad)".
- Right, "MEASURED SLOPE": a = 1.84 ± 0.02 rad/(rad/s), regression and calibration line.
- "PREDICTED 8π²NR²/(λc)": a_th = 1.7 ± 0.1.
- "AGREEMENT": a and a_th with their error bars on one axis. "0.7σ apart, p = 0.47, with the regression error only".
- "LINEAR REGRESSION": ten setting means, equal weights, R² = 0.9986, ω² term: p = 0.65, no bend.

**Say.** "Each point is the mean phase of one setting, against its calibrated ω. The points lie on a straight line, and the differences below show no pattern. A linear regression gives a slope of 1.84 ± 0.02 radians per rad/s. The prediction from the given loop parameters is 1.7 ± 0.1. On the right both are drawn with their error bars: they overlap, and the difference is 0.7 standard deviations, p = 0.47. The prediction's error bar is much wider than ours, because of the turn count."

**In depth.**

*The fit.* $\Delta\phi = a\omega + b$, an ordinary linear regression through the ten setting means, all with equal weight. This is how Schwartz obtains his slopes: "fitting a linear regression model to the data (Δφ = a·Ω + b)", quoting the regression's standard errors [Schwartz Sec. 3, Eq. 7].

*Where ± 0.02 comes from* (Part 1.5):
- The regression measures how far the points lie from the line and uses it as the uncertainty of every point: $s = \sqrt{\sum(\Delta\phi - \text{line})^2/(10 - 2)} = 0.045$ rad.
- The slope's standard error is $u_a = s/\sqrt{\sum(\omega - \bar\omega)^2} = 0.024$.
- The calibration line's own uncertainty (k and ω₀ are not exact) adds 0.004: it is propagated by redrawing the line's two coefficients 4000 times and refitting. In quadrature: 0.024 again, shown as ± 0.02.

*R² = 0.9986* is the fraction of the variation of Δφ that the line accounts for. In a regression χ²/dof = 1 by construction, so R² is the fit-quality number here (on slide 4 the points had error bars of their own, so χ² was used there).

*Why 1.84 and 1.7 agree although they look far apart.* Both are rounded. Unrounded, a = 1.839 ± 0.024 and a_th = 1.738 ± 0.137. The difference is 0.101, not 0.14. Combined uncertainty $\sqrt{0.024^2 + 0.137^2} = 0.139$, so $z = 0.101/0.139 = 0.72$, shown as 0.7σ. Two-sided $p = 2(1 - \Phi(0.72)) = 0.47$: the probability of a difference at least this large if the two were equal. "With the regression error only": slide 9 adds the method and calibration-model uncertainties, which give 0.7σ and p = 0.51 (z = 0.66).

*Why not use the photos' own error bars.* Each setting mean has a standard error from its photos (SEM, 0.004–0.028 rad). Weighting by SEM and the table-speed term gives χ²/dof = 56.6/8, p = 2 × 10⁻⁹: the points lie 2.7 times farther from the line than those error bars allow. The SEM assumes the photos of one setting are independent, and they are not: the phase wanders slowly during the run (B4). The regression's s contains every effect that moves a setting mean (photo noise, table speed, wander), so it is the honest per-point uncertainty. Our written analysis instead added an extra per-setting term σ_x = 0.041 rad chosen so that χ²/dof = 1, and got the same 1.839 ± 0.024 [`sagnac_update/sagnac.txt`].

*Two checks of where the extra variation comes from* (computed for the presentation):
1. Split each setting's photos into first and second halves. The halves differ by more than their photo noise allows in 8 of 10 settings (2.5–6 standard errors). From these differences a setting mean is expected to wander by 0.039(7) rad, close to s = 0.045 rad.
2. The fringe spacing P (mean per setting) changes by up to ±2 % between settings (ring fit and cosine fit agree, r = 0.89), and the points' distances from the line follow it, r = +0.81, p = 0.005. Something changes the two beams slightly between settings, and the phase moves with it (the fiber guides several modes, Part 1.4).

**Numbers and where they come from.** Setting means: `sagnac_update/sagnac_speeds.csv` (also B3). Fit: `slides.py` `fit8('regression')`. Prediction: slide 3.

**Questions you may get.**
- *Why no error bars on the points?* The regression gives every point the same uncertainty, s = 0.045 rad, which is smaller than the markers.
- *Why does the line miss the origin?* The intercept is −0.05 ± 0.03 rad, 1.8 standard errors from zero (p = 0.11), not significant. Sagnac predicts zero, but the rest point itself is uncertain: 10 photos at the start, and setting means wander by about 0.045 rad. Forcing the line through the origin gives a = 1.80 ± 0.01. That 0.04 difference is smaller than the method uncertainty of slide 9, and both values agree with the prediction.
- *Is it really a straight line?* That is the last line on the slide: add an ω² term to the fit, and it comes out 0.5 standard errors from zero (p = 0.65, t test). No sign of a bend. R² alone would not show this, since almost any smoothly rising data gives R² near 1.
- *Why R² here when slide 4 insists on χ²?* Both slides first tried χ² with the points' own error bars, and both failed (47.2/7 on slide 4, 56.6/8 here). The difference is what made them fail. On slide 4 the points bend in a pattern that a quadratic removes, a real shape, so we keep the χ² result and carry the bend as an uncertainty. Here the points jump up and down with no pattern (ω² term p = 0.65): extra random variation between settings, not a wrong model. Then a regression, which takes the noise from the points themselves, is the standard tool (Schwartz Eq. 7). In a regression χ²/dof = 1 by construction, so R² describes the fit and the ω² term tests the shape.
- *Why not fit all 161 photos instead of ten means?* That would treat correlated photos as independent and give 1.85 ± 0.01, more than twice as precise as the regression's ± 0.02, which it is not.
- *How does this compare with Schwartz?* His large ring gave 1.304 ± 0.038 against a prediction of 1.344 ± 0.006, agreement within 3 % (his three cases agree within 3–8 %). Our a/a_th = 1.06 ± 0.08, with a much less precise prediction because of N.
- *What is the rest point's error?* Its SEM is 0.02 rad, and it is common to every other point, a shared offset that the free intercept absorbs.

**Transition.** "That slope used one way of reading the phase from the photos. Does the answer depend on the method?"

## Slide 09. Three full-image methods differ by at most 6 %

**Purpose.** Test the phase extraction with other methods, and build the final uncertainty of the slope.

**On the slide.**
- Subtitle: the same 161 photographs analyzed with four ways of reading the phase.
- Left, "Phase change relative to the ring-minima result": for each setting, the cosine fit's and template matching's phase change minus the ring minima's.
- Right, "Fitted slope (gray: prediction ± 1σ)": the slopes 1.84 (ring minima, main), 1.88 (cosine fit of the profile), 1.77 (template match to rest profile), with their error bars and a bracket marking their range. Below a dashed line, "10 random pixels (mean of 1000)": 1.70. The gray band is the prediction ± 1σ.
- Bold: a = 1.84 ± 0.07 rad/(rad/s). Under it: 1.3 % regression and calibration line, 3.0 % method, 1.7 % calibration model, in quadrature 3.7 %.
- Footer: each slope comes from the same linear regression as the previous slide, R² = 0.9986, 0.9982, 0.9993 from top. Error bars: regression and calibration line.

**Say.** "Does the result depend on how we read the phase? We applied two more methods to the same photographs: a fit of a cosine to each photo's profile, and matching each profile to the rest profile. The three methods that use the whole image differ by at most 6 %, and we take that spread as a method uncertainty. Together with the regression and the calibration this gives our final slope, 1.84 ± 0.07. The fourth row, ten random pixels, is the method on the next slide."

**In depth.**

*The methods.*
- **Ring minima** (main, slides 6–7): the center phase from the ring radii. Needs no model of the brightness.
- **Cosine fit** [Prep Q5 asks for two ways to extract Δφ]. The intensity against $s = r^2$ is fitted with $C e^{-s/s_C} + B e^{-s/s_B}\cos(2\pi s/P + \phi) + D$. The exponentials are the envelopes of two Gaussian beams of different size [QE3, PDF p. 7], D is the dark level, φ the center phase. It uses every pixel, but relies on the envelope model and a linear camera response.
- **Template matching** (our own method, not from the documents). In words: slide the at-rest profile along r² until it lies on top of the new one, and read the phase from how far you had to slide it.
  1. *The template.* Each photo is turned into a profile: brightness against $s = r^2$ (averaged over the fringe sector, as on the right of slide 6). The 10 rest profiles are averaged into one template $T(s)$.
  2. *The match.* For each photo, find the shift $D$ (in px²) for which $g\,T(s + D) + c$ best fits that photo's profile over $s$ = 0.2–1.8 × 10⁵ px² (least squares). The gain $g$ and offset $c$ absorb changes of overall brightness and background.
  3. *Shift to phase.* Shifting by one fringe spacing $P$ in $s$ is one full fringe, so $\Delta\phi = 2\pi D/P_T$, where $P_T$ is the template's own spacing (the shift that maps the template onto itself one fringe further out). Rings that moved inward give $D > 0$.
  - Strength: no model of the fringe shape, it uses the whole measured shape. Weakness: one shift describes the whole window, so it reads an average phase over the window, effectively near its middle, not at the center. If the ring spacing $P$ changes between photos, the rings far out move more than the center ones, and that adds to the phase read mid-window, but not to the center phase. This is why its per-setting phases differ from the minima by up to 0.16 rad.

*Left panel.* The cosine fit reads 0.10–0.19 rad higher than the minima at every moving setting. The difference comes from the rest photos: the bright, flat center at rest is not a cosine (fit rms 9.5 counts at rest, 1.4–2.3 counts above S = 44) [`sagnac_update/sagnac.txt`]. The cosine fits report an rms and no χ², because neighboring pixels share speckle, so the profile points are not independent. Template matching differs from the minima by −0.16 to +0.07 rad from setting to setting.

*Method uncertainty.* The three methods read the same photographs, so they are correlated checks, not independent measurements. Their differences show how much the answer depends on the choice of method. Half the range, (1.876 − 1.768)/2 = 0.054 rad/(rad/s), 3.0 % of a, is taken as one standard uncertainty. "At most 6 %": the slopes span 0.11, 6 % of 1.84. Their error bars do not all overlap, which is exactly why the difference is counted as an uncertainty and not ignored.

*The final slope.* Three parts, added in quadrature:

| Part | On a | Relative |
|---|---|---|
| Regression and calibration line (slide 8) | 0.024 | 1.3 % |
| Phase-reading method | 0.054 | 3.0 % |
| Calibration model, linear vs quadratic (slide 4, B1) | 0.032 | 1.7 % |
| Total | 0.067 | 3.7 % |

So a = 1.84 ± 0.07 rad/(rad/s). Against the prediction: 0.7σ, p = 0.51.

*Why quadrature, and how it could be attacked.* Independent sources add in quadrature because their variances add. The three parts come from different things: noise between settings, the choice of phase reading, the choice of calibration curve. They are close to independent: the methods share the same photos and calibration, so the regression noise is mostly common to them and drops out of their differences, and the calibration model shifts ω for all methods alike, so it does not enter the method range.
- *"Half the range as one σ is arbitrary."* It is a convention. Range/√12 (a flat distribution over the range) would give 0.031, the full range 0.108. With the full range the total is 0.115 and the difference from theory 0.6σ.
- *"Why not add them linearly?"* That assumes all three push the same way at full size, the worst case: 0.024 + 0.054 + 0.032 = 0.110, again 0.6σ. No choice changes the conclusion, because the prediction's own ± 0.14 is the largest uncertainty.
- *"Is the regression error counted twice?"* A little regression noise may also sit in the method differences. If so, the total is slightly too large, the cautious direction.

*The pixel row.* This is the instructor's suggested method (slide 10). It shows the mean over 1000 random draws, 1.70, not one draw. Its error bar is 0.06, the typical regression error of one draw. The mean's own uncertainty (0.006) would be too small to show here, because all draws use the same 161 photographs, so the wander of the setting means is common to all draws and does not average out. The pixel method is not used for the result or the method uncertainty (slide 10).

**Numbers and where they come from.** Per-setting phases of each method: `sagnac_update/sagnac_speeds.csv`. Slopes refitted with the same regression in `slides.py` (so the last digit can differ from `sagnac_update/RESULTS.md`, which used a weighted fit). Linear vs quadratic: B1.

**Questions you may get.**
- *Which method is right?* None is "right" by itself. The minima method is the main one because it reads the phase at the center and needs no brightness model. The spread between methods is the uncertainty.
- *Do the methods also fail the photo-error test?* Yes. With each method's own photo errors (SEM) and the table speed, χ²/dof = 56.6/8, 131.6/8 and 42.7/8, p = 2 × 10⁻⁹, 1 × 10⁻²⁴ and 1 × 10⁻⁶. All three see the same extra variation between settings, so it is in the photos, not in one method. Each method's points lie typically 0.045, 0.053 and 0.031 rad from its line.
- *Your written report says ± 0.06.* It lists the calibration-model change separately instead of adding it. The deck adds it, hence ± 0.07.
- *Where is the 2D Fourier method of your earlier report?* The final analysis keeps three methods. (Have your reason ready if asked. The earlier report gave 1.90 for it with the quadratic calibration. Refitted with the linear calibration and the same regression it gives 1.86, inside the range of the three, so keeping it would not change the method uncertainty.)

**Transition** (if showing slide 10). "The fourth row is the method with random pixels. Here is what happens when we repeat it."
**Transition** (if skipping it). "Now, what does this slope tell us about the speed of light?"

## Slide 10. The slope from ten random pixels depends on the draw (optional)

**Purpose.** Show how much the random-pixel method depends on the pixels drawn, fairly and with its reasons. This is the method the instructor suggested, so present it as a test of the method under our conditions, not as a criticism: the main reason it struggles is a property of our run (half a fringe).

**On the slide.**
- Subtitle: the pixel method repeated for 1000 random draws of ten pixels.
- Left: histogram of the slope over the draws (gray: fixed pixel positions, the method itself). Purple outline: the same pixels, but each photo's own fitted center and P used. Vertical lines: ring minima 1.84 (blue), prediction 1.7 ± 0.1 (gray dashed).
- Right, "SPREAD OVER 1000 DRAWS": mean slope 1.702(6), spread (SD) 0.196(5), typical fit error 0.057, spread ÷ fit error 3.4.
- "MORE PIXELS PER DRAW": the spread between draws and one draw's fit error for 10, 20, 40, 80 pixels.
- "LIMITATIONS IN OUR RUN": the rings move by only half a fringe, so each pixel sees less than one bright-dark cycle. A fixed pixel also reads a moving center or a change of P as a phase shift.
- Footer: each draw is 10 pixels in the fringe sector, read at the same image position in all 161 photos. One least-squares fit gives every photo's phase. Spread: standard deviation of the slopes. 2 draws with a negative slope (failed fits) are left out.

**Say.** "This is the method of following individual pixels: pick ten pixels in the fringe region and follow their brightness through all photos. As the rings move, each pixel gets brighter and darker, and one fit finds the phase of every photo from all ten pixels together. The question is how much the answer depends on which ten pixels you pick. We repeated the whole method a thousand times. One draw typically gives a slope between about 1.5 and 1.9: a spread of 0.2, three to four times larger than the error bar a single draw reports. The main reason is our run: the rings move only half a fringe, so each pixel sees less than one bright-dark cycle."

**In depth.**

*The model* [Analysis, `sagnac_update/methods.py` header]. For pixel j in photo k:
$$I_{jk} = a_j + b_j\cos(\psi_j + \delta_k)$$
$a_j$ is an offset, $b_j$ an amplitude, $\psi_j = 2\pi s_j/P_0$ the pixel's fringe phase at rest ($s_j$ its r² from the rest ring center, $P_0$ the rest spacing), and $\delta_k$ the photo's phase change. An alternating least-squares fit (Part 1.5: it fits the pixel parameters with the phases fixed, then the phases with the pixel parameters fixed, and repeats) finds all $a_j$, $b_j$ and $\delta_k$, with the physical limits $b_j \ge 0$ and $a_j - b_j \ge$ dark level. It is the same two-beam interference law as Manual Eq. 6. Each draw's phases go through the same averaging per setting and the same regression as slide 8.

*The scan.* 1000 draws of ten pixels (seeds 1–1000), each exactly like the single draw of our written analysis (seed 20261011, slope 1.79, well inside the histogram) [`presentation_slides/source/random_points_scan.py`].

*The table.*
- Mean slope 1.702(6): the average of the 998 good draws. The (6) is SD/√998, how well this average is known.
- Spread (SD) 0.196(5): how much one draw differs from another. The (5) is from resampling the draws (bootstrap). A robust spread (half the width of the middle 68 %, 1.52–1.89) gives 0.187, almost the same.
- Typical fit error 0.057: the median of the slope errors each draw's own regression reports.
- Spread ÷ fit error 3.4: one draw's error bar understates its real uncertainty by a factor of about 3.4. One draw of ten pixels gives the slope only to about ±0.2, about 11 %.

*Failures.* 2 of 1000 draws gave a negative slope: failed fits, left out. 18 draws stopped at the 500-round limit of the fit, and their slopes are kept.

*Limitations in our run.*
1. **Half a fringe.** Our largest shift is 0.52 fringe, so each pixel passes through less than one bright-dark cycle. The fit must tell a phase change from a change of brightness or contrast with little to go on, which is why one draw depends so much on its pixels. This is a property of our run, not of the method. Prep Q5 asks which way of reading the phase is more accurate when Δφ_max is large: with a shift of several fringes each pixel passes through full cycles and following fixed points works well. Our shift is small, the case where ring positions do better.
2. **Fixed positions.** A pixel at a fixed position reads a moving ring center or a change of P as a phase change. The purple diagnostic computes each pixel's phase from the photo's own fitted center and spacing: the mean moves from 1.70 to 1.78 and the spread narrows to 0.17. So about half of the 0.14 difference from the ring minima comes from the fixed positions.

*More pixels.* Spread 0.196, 0.091, 0.073, 0.056 for 10, 20, 40, 80 pixels (60 draws each above 10), while one draw's fit error is 0.057, 0.041, 0.040, 0.037. The spread shrinks but stays above the fit error, and the mean settles near 1.73 (1.70, 1.71, 1.73, 1.73). More pixels reduce the luck of the draw, not the effect of fixed positions.

**Questions you may get.**
- *"The pixel method gives 1.70, closer to the prediction than your 1.84."* The prediction's ± 0.14 cannot tell 1.70 from 1.84: both are within one standard deviation of it. The point of the slide is not which mean is closer, but that one draw of ten pixels gives anything from about 1.5 to 1.9, while its own error bar says ± 0.06.
- *Would the method work better in another run?* Yes, with a larger phase shift (more turns, faster rotation), each pixel would see full cycles. And with more pixels.
- *Why is the histogram not centered on 1.84?* Partly fixed positions (purple), partly the small shift.

**Transition.** "Back to the main result: what does the slope tell us about the speed of light?"

## Slide 11. The slope gives the speed of light, 0.7σ from c₀

**Purpose.** Turn the slope into c, show the result is consistent with the defined value, and show what limits it.

**On the slide.**
- Subtitle: the Sagnac slope solved for c, with the measured a and the given N, R and λ.
- Formula and result: $c = 8\pi^2 N R^2/(\lambda a) = (2.8 \pm 0.2) \times 10^8$ m/s.
- "COMPARED WITH THE DEFINED VALUE": our c with its error bar on a number line, and c₀ dashed. "0.7σ from c₀, p = 0.50".
- "UNCERTAINTY BUDGET (ADDED IN QUADRATURE)": turn count 7.7 %, phase-reading method 3.0 %, loop radius 1.8 %, calibration model 1.7 %, regression and calibration line 1.3 %. Total of the listed sources 8.7 %.
- "Not included: a steady drift of the phase in time, and λ (given without an uncertainty)."
- Footer: c₀ = 299 792 458 m/s (defined).

**Say.** "The predicted slope contains the speed of light. If we trust N, R and λ, we can turn our measured slope around and solve for c. We get 2.8 ± 0.2 times 10⁸ meters per second, 0.7 standard deviations from the defined value. The uncertainty is 8.7 %, and most of it is the given turn count, 13 ± 1. Everything from our own measurement together is 3.7 %. Not included is a possible steady drift of the phase in time, which I come back to on the next slide."

**In depth.**
- From $a = 8\pi^2 N R^2/(\lambda c)$ [Manual Eq. 5 with A = πR²]: $c = 8\pi^2 N R^2/(\lambda a)$.
- Numbers: $c = 8\pi^2 \cdot 13 \cdot (0.570\ \mathrm{m})^2/(640\ \mathrm{nm} \cdot 1.839\ \mathrm{s}) = 2.834 \times 10^8$ m/s, uncertainty 0.246 × 10⁸, rounded (2.8 ± 0.2) × 10⁸ m/s. (rad/(rad/s) is a time, so a is in seconds.)
- c/c₀ = 0.95 ± 0.08, $z = -0.67$, two-sided p = 0.50.
- **Budget.** $c \propto N R^2/a$, so the relative uncertainties add in quadrature, R counted twice:
  $$\left(\frac{u_c}{c}\right)^2 = \left(\frac{u_N}{N}\right)^2 + \left(2\frac{u_R}{R}\right)^2 + \left(\frac{u_a}{a}\right)^2$$
  with $u_a/a$ split into its three parts (slide 9). N: 1/13 = 7.7 %. Method 3.0 %. R: 2 × 0.88 % = 1.8 %. Calibration model 1.7 %. Regression and calibration line 1.3 % (regression 1.3 %, calibration line 0.2 %). Total 8.7 %.
- **Not included.** A steady drift (slide 12, B4). The wavelength: λ was given without an uncertainty, and each nanometer would change c by 0.16 %.
- **Not an independent test.** $c/c_0 = a_{\mathrm{th}}/a$: this is the same comparison as slide 8, shown as a speed of light. And c is inferred through the Sagnac formula, not timed.

**Questions you may get.**
- *Why is c lower than c₀?* It is 0.7σ low, well within the uncertainty. The same 0.7σ as the slope being high.
- *Turned around, what does the slope say about N?* With c, R and λ taken as known, N = 13.8 ± 0.5 turns (± 0.6 including R), about twice as precise as the given 13 ± 1. Counting the turns would be the single best improvement.
- *What if 13 ± 1 is a bound?* Then u(N) = 0.58 turns and the total drops to about 6 %.
- *Could you do better with the laser frequency?* The phase fixes the product λc. With a known frequency ν instead of λ, $c = \sqrt{8\pi^2 N R^2 \nu/a}$, with half the relative uncertainty.
- *Why does slide 8 say p = 0.47 and this slide p = 0.50?* Slide 8 uses only the regression error of a. Here the method and calibration-model parts are included as well. With those parts the slope comparison itself gives z = 0.66, p = 0.51. For c the same uncertainties are propagated as relative uncertainties (c ∝ 1/a, to first order), which gives z = −0.67, p = 0.50. The two differ only in the second digit.

**Transition.** "To sum up."

## Slide 12. The Sagnac phase is linear in ω and matches theory

**Purpose.** The three results, the main limitation, and what to do next.

**On the slide.**
- 1.84 ± 0.07 rad/(rad/s): slope of Δφ against ω, full uncertainty. Prediction 8π²NR²/(λc) = 1.7 ± 0.1, 0.7σ apart.
- 3.29 ± 0.06 rad (0.52 fringe): largest shift, at ω = 1.81 rad/s. A travel-time difference of only 1.1 fs.
- (2.8 ± 0.2) × 10⁸ m/s: inferred c, 0.7σ from the defined value. Limited by the turn count N = 13 ± 1.
- "LIMITATION": a slow drift of the phase in time would change the slope, since the speeds only went up.
- "NEXT": rest photos between the speeds and reversed rotation, count the fiber turns, couple the light into a single fiber mode (possible with this fiber, but hard to align).

**Say.** "The phase grows linearly with the rate of rotation, with a slope of 1.84 ± 0.07 radians per rad/s, 0.7 standard deviations from the prediction. At our top speed that is 3.29 radians, half a fringe, a travel-time difference of about one femtosecond. Read as a speed of light, 2.8 ± 0.2 times 10⁸ meters per second, limited by the turn count. Our main limitation: we took the speeds in increasing order, so a slow drift of the phase in time would change the slope and we cannot separate it. Next time we would take rest photos between the speeds and reverse the rotation, count the turns, and couple the light into a single fiber mode."

**In depth.**
- **3.29 ± 0.06 rad.** Δφ at S = 90 is the difference of two setting means (top speed minus rest), each uncertain by about s = 0.045 rad, so $\sqrt2 \times 0.045 = 0.06$. In time, $\Delta t = \Delta\phi\,\lambda/(2\pi c) = 1.1$ fs, inferred from the phase, not timed. Predicted: 3.15 rad, 1.07 fs.
- **The limitation.** Rest photos only at the start, speeds only going up: ω and time rise together, so a phase drifting steadily in time looks like extra slope (or less, depending on its sign). The phase does drift within single settings, by −6.9 to +8.2 mrad/s, but the sign changes from setting to setting. Their plain mean, −1 ± 1 mrad/s, cannot rule a steady drift in or out. (A mean weighted by their own error bars would be misleading, because the ten drifts differ far more than those error bars allow.) A joint fit cannot separate the two either (B4). So the quoted uncertainty does not cover a steady drift.
- **Next.**
  - Rest photos between the speeds and reversed rotation: the Sagnac phase changes sign with the rotation direction, a drift does not. Sagnac himself measured this way, in one direction and then the other [Pascoli Sec. 2]. Random order of speeds also helps.
  - Count the turns: N = 13 ± 1 is the largest uncertainty.
  - Single-mode coupling: at 640 nm the fiber guides several modes (V = 5.6). The manual computes up to 93 % in the fundamental mode with a matched 15.5 mm lens and at most 47 % with the setup's 8 mm lens [Manual p. 10]. It is possible but was very hard to align. With mostly one mode the light should be less sensitive to stress and temperature, a plausible source of the wander (Part 1.4).

**Questions you may get.**
- *Did you check the sign of the effect?* No, the rotation direction was not recorded.
- *How sure are you the wander comes from the fiber modes?* Not sure: it is a plausible explanation, not tested. What we know is that the phase wanders more than the photo noise and follows the fringe spacing (slide 8).

**Transition.** "Thank you."

## Slide 13. Thank you

**On the slide.** "Thank you", "Questions?", names, instructor, and one photo at ω = 0.88 rad/s.

Open for questions. Use the backup slides (B1–B4) when a question matches them. For questions about the result, go back to slide 12.

---

## Backup B1. A quadratic calibration would raise the slope by 1.7 %

**When to open it.** "Why a straight line when the calibration points bend?" or "How much does the calibration model matter?"

**On the slide.**
- Left: measured ω minus each model, in % of ω, against S, for the linear fit (used, orange) and the quadratic (purple). The linear residuals form an arch, falling to −2.3 % at S = 90. The quadratic ones lie around zero within their error bars.
- Right, a table, linear vs quadratic: χ²/dof 47.2/7 vs 6.0/6, p 5 × 10⁻⁸ vs 0.42, slope a 1.839 vs 1.871, a/a_th 1.058 vs 1.076, c 2.834 vs 2.786 (10⁸ m/s), ω at S = 90 1.814 vs 1.782 rad/s.
- Text: the main analysis uses the linear model, but the table turns slower than linear at high settings. The choice moves a and c by 1.7 %, included as the calibration-model uncertainty.
- Footer: both fits include ε = 0.48 %, set by χ²/dof = 1 for the quadratic, so its p is not an independent test. A cubic adds nothing (χ²/dof = 2.6/5, p = 0.75).

**Say.** "The straight line misses the top points by up to 2.3 %. A quadratic follows them. Using it would raise our slope and lower c by 1.7 %. Since we cannot be sure which calibration is right for the photos, that 1.7 % is in our uncertainty."

**In depth.**
- The task asks for the linear relation [Manual minimum task 1]. The quadratic's χ²/dof = 6.0/6 is by construction: ε was chosen to make it so. Its p is therefore not evidence for the quadratic. The cubic adds nothing.
- Same regression as slide 8 with the quadratic ω: 1.871, +1.7 %. c: 2.786 × 10⁸ m/s, −1.7 %.
- Calibration-model term: |1.871 − 1.839| = 0.032 on a, 1.7 %.
- **If the instructor has seen our earlier report.** The first version of our report used this quadratic calibration as its main result: a = 1.87 ± 0.06 with the method term, 0.9σ from the prediction, c/c₀ = 0.93 ± 0.08, largest shift 3.29 ± 0.07 rad at ω = 1.78 rad/s. It also listed a 2D Fourier method (1.90 then, 1.86 with the linear calibration) and a phase-shifting method, AIA (1.61, not used). The final analysis uses the linear calibration, keeps three phase-reading methods, and carries the linear-to-quadratic difference as an uncertainty. The old report is in `previous_reports/`.

## Backup B2. Uncertainty sources considered

**When to open it.** "Which uncertainties did you include?" This is manual task 5 in one table.

**On the slide.** A table with source, kind, size and whether it is included:

| Source | Kind | Size | Included |
|---|---|---|---|
| Photo-to-photo variation | statistical | SEM 0.004–0.028 rad | yes, regression |
| Table speed at a setting | statistical (calibration) | ε = 0.48 % of ω | yes, regression |
| Wander between settings | statistical | about 0.04 rad | yes, regression |
| Calibration line | calibration, correlated | 0.2 % of a | yes |
| Phase-reading method | systematic | 3.0 % of a | yes |
| Turn count N = 13 ± 1 | given | 7.7 % | yes (theory, c) |
| Loop radius R = 57.0 ± 0.5 cm | given | 1.8 % | yes (theory, c) |
| Wavelength λ = 640 nm | not given | 0.16 % per nm | no |
| Steady drift in time | systematic | not separable from a | no |
| Linear vs quadratic calibration | model choice | 1.7 % of a | yes |

Subtitle: N, R and λ enter the prediction and c. All other sources enter the measured slope a. Footer: SEM defined, statistical parts enter through the regression, others in quadrature, a = 1.84 ± 0.07. Not checked: the sign of the shift. ω comes from the calibration, not measured during the photos.

**In depth.**
- "Yes, regression": photo noise, table speed and wander all move the setting means away from the line, and the regression takes its per-point uncertainty from exactly that distance (slide 8). So they are in the ± 0.024 without being added separately.
- A steady drift would **not** show up as distance from the line, because it would move all points along with ω (it changes the slope instead). That is why it is the one effect the regression cannot see.
- The calibration line, the method and the calibration model are added in quadrature afterwards (slide 9). N and R enter only the prediction and c.

## Backup B3. Measured values at each setting

**When to open it.** "Can we see the numbers behind slide 8?"

**On the slide.** For each setting: motor setting S, number of photos, time of the first photo, ω (from the linear calibration), Δφ (mean, with SEM) and the ring spacing P:

| S | Photos | First photo | ω (rad/s) | Δφ (rad) | P (10⁵ px²) |
|---|---|---|---|---|---|
| rest | 10 | 16:21:40 | 0 | 0 (reference) | 1.075 |
| 12.59 | 15 | 16:24:10 | 0.235(1) | 0.39(3) | 1.071 |
| 20.74 | 15 | 16:25:44 | 0.402(2) | 0.64(2) | 1.058 |
| 31.11 | 19 | 16:26:38 | 0.613(3) | 1.02(2) | 1.065 |
| 44.07 | 22 | 16:27:54 | 0.877(5) | 1.59(2) | 1.093 |
| 55.83 | 18 | 16:29:08 | 1.117(6) | 1.984(8) | 1.058 |
| 64.00 | 16 | 16:30:16 | 1.283(7) | 2.375(8) | 1.094 |
| 73.65 | 16 | 16:31:28 | 1.480(8) | 2.700(4) | 1.087 |
| 82.18 | 15 | 16:32:32 | 1.654(9) | 2.958(4) | 1.068 |
| 90.00 | 15 | 16:33:24 | 1.814(9) | 3.285(5) | 1.073 |

**In depth.**
- 0.235(1) means 0.235 ± 0.001.
- The Δφ errors are the photo-to-photo SEM. The regression does not use them (slide 8): the setting means vary more than these allow.
- The rest row is the reference, 0 by construction. Its SEM, 0.02 rad, is shared by every other row: a common offset, not independent noise on each point.
- ω errors combine the calibration line and the 0.48 % reproducibility.
- **3.285(5) here and 3.29 ± 0.06 on slide 12.** The (5) is the SEM of the top setting's photos only. Slide 12's ± 0.06 is the realistic uncertainty of a difference of two setting means, each uncertain by s = 0.045 rad.

## Backup B4. A steady phase drift would mimic part of the slope

**When to open it.** "Could something other than rotation change the phase during the run?" or any question about drift or the order of measurements.

**On the slide.**
- Top: ω against minutes since the first photograph, a staircase that only goes up.
- Bottom: each photo's phase minus the fitted line of slide 8, against time, with a straight line through each setting's photos.
- Right, "JOINT FIT": $\Delta\phi = a\omega + Dt + b$, a = 2.1 ± 0.2 rad/(rad/s), D = −0.8 ± 0.5 mrad/s, correlation(a, D) = −0.99, drift term p = 0.12, R² = 0.9991. Without D: a = 1.84 ± 0.02.
- "WITHIN EACH SETTING": trend −6.9 to +8.2 mrad/s, 6 of 10 settings with p < 0.05.
- Footer: at 640 nm SMF-28e+ guides several modes (V number 5.6, above the single-mode limit 2.405). Unless the light is coupled into one mode, it travels in several, and light in several modes is sensitive to stress and temperature, a plausible source of the wander. The quoted uncertainties do not include a steady drift.

**Say.** "We photographed the speeds in increasing order, with rest photos only at the start, so ω and time rose together. A phase that drifts steadily in time would look exactly like extra slope. Fitting both a slope and a drift gives a drift consistent with zero, but the two fitted values have a correlation of −0.99: the data cannot tell them apart. Within single settings the phase does drift, by a few mrad per second, but the sign changes from setting to setting. So we report the slope without a drift term and say clearly that a steady drift is not covered."

**In depth.**
- **Joint regression** with t the mean time of each setting's photos, as in slide 8. D is not significant (t test, p = 0.12). The correlation −0.99 means raising a and lowering D fits almost equally well. If D were real, being negative, it would have lowered the slope, so the slope without D would be too small, not too large. With the first photo's time instead of the mean, D = −1.1 ± 0.5, p = 0.07, the same conclusion.
- **Within settings** the drifts are real (6 of 10 significant) but change sign. Their plain mean, −1 ± 1 mrad/s, cannot pin down a steady drift. Between settings the means wander by about 0.045 rad, the s of slide 8.
- **A plausible source** (not tested): the fiber guides several modes at 640 nm (Part 1.4). Light in several modes is sensitive to fiber stress and temperature, and the camera sees the whole multimode field, so the two directions need not stay exactly reciprocal.
- **How to separate it next time.** Rest photos between the speeds, and reverse the rotation. The Sagnac phase changes sign with the rotation direction, a drift does not [Pascoli Sec. 2]. Random order of speeds.

**Questions you may get.**
- *Both beams go through the same fiber, so why does stress not cancel?* For a single mode it does (reciprocity). With several modes the power exchanged between modes along the fiber is not the same for the two directions, and the camera sees the whole multimode field, so the cancellation is not guaranteed. Fiber gyroscopes put a single-mode filter at the common port for exactly this reason (general knowledge, not in our documents).
- *So is your slope wrong?* We cannot tell. The joint fit says a drift is not needed (p = 0.12), and it says the data cannot separate the two either. That is why it is stated as a limitation and not hidden in the uncertainty.

---

# Part 3. Question bank

Short answers, sorted by topic. The slide in brackets is where the full answer is.

**Physics**

- *Why does rotation delay one beam?* The splitter moves while the light travels, so the co-rotating beam has farther to go (slide 2).
- *Why does a delay show up as a phase?* Both beams carry the same light wave. Delaying one by Δt shifts it by Δt/T of a period, a phase 2πΔt/T, and that moves the interference bands (slide 2).
- *Why do the rings move inward rather than just change brightness?* Rotation adds the same Δφ everywhere. Each dark ring must stay at φ = π, 3π, …, so its r² drops by PΔφ/2π, the same for every ring (say it on slide 6).
- *Is this a relativistic effect?* It is first order in v/c, where the classical and relativistic results agree. It needs only that light moves at c both ways in the lab frame (Einstein's second postulate). Corrections are of order 10⁻¹⁷ (Part 1.1).
- *Does the refractive index of the fiber matter?* No. With Galilean velocity addition a fiber loop would show no effect at all. Relativistic addition (Fresnel drag) removes n and gives the vacuum result (Part 1.1) [Schwartz Sec. 2, Prep Q10].
- *Does the position of the rotation axis matter?* No, only the enclosed area (Part 1.1) [Schwartz Eq. 4, Prep Q8].
- *Why rings and not straight fringes?* The beams are Gaussian beams with different wavefront radii (slide 6, Part 1.3).
- *Why is the effect so small?* $\Delta t = 4NA\omega/c^2 \approx 1$ fs. It is measurable because interference compares the beams to a fraction of the 2.13 fs optical period (slide 2).
- *Which way should the rings move?* That depends on the rotation direction and on which beam has the larger wavefront radius. We did not record the direction, so we compare magnitudes only (slide 12).

**Setup and calibration**

- *Where do N, R and λ come from?* Given in the manual's setup figure. We did not count the turns (slide 3).
- *What do the polarizers do?* Suppress reflections that reach the camera without going around the loop (slide 3, Part 1.2) [Prep Q6].
- *Is the fiber single-mode?* Not at 640 nm: V = 5.6 > 2.405 (B4, Part 1.4) [Prep Q12].
- *SMF-28e or SMF-28e+?* The manual's text says SMF-28e, its printed datasheet and the preparation sheet say SMF-28e+. Same core and NA as used here (slide 3).
- *What is S?* The number set on the motor controller, which sets the motor voltage. No unit (slide 4).
- *Why a straight-line calibration if it fails χ²?* The task asks for it, and the line is good to 2.3 %. It fails only because the error bars are 0.5 %. The 1.7 % a quadratic would change is in the uncertainty (slide 4, B1).
- *How is u(T) from stopwatch splits derived?* Only the first and last press set the mean turn time (slide 4).
- *Was ω measured during the photos?* No, from the calibration at the same setting. ε covers the variation (slide 4).
- *Why was task 4 (the small off-center loop) not done?* It was skipped on the instructor's instruction.

**Phase extraction**

- *How do you get one phase per photo?* Fit $R_m^2$ against ring number: the intercept gives η, $\phi_c = \pi - 2\pi\eta$ (slide 7).
- *Do you need the pixel size?* No, η is a ratio of two px² numbers (slide 7).
- *How good are the circle fits?* Single minima within 3–8 px (rms), radii to about ±1 px (slide 6).
- *Why the center phase?* Rotation changes only φ_c, and a change of ring spacing does not move it (slides 6–7).
- *Why does the cosine fit read higher?* Its model does not describe the bright flat center at rest (slide 9).

**Fit and statistics**

- *Why a regression and not a weighted fit?* The setting means vary more than the photo errors allow (χ²/dof = 56.6/8). The regression uses the actual distance of the points from the line (slide 8).
- *Why R² on slide 8 but χ² on slide 4?* Both failed χ² with their own error bars. Slide 4's points bend in a pattern (a real shape, carried as an uncertainty). Slide 8's vary randomly with no bend (ω² term p = 0.65), so a regression that takes the noise from the points is right, and there χ²/dof = 1 by construction (slide 8, Part 1.5).
- *What does p = 0.47 mean?* The chance of a difference at least this large if the slope equaled the prediction. Not the chance that the theory is right (Part 1.5).
- *Why 0.7σ if 1.84 and 1.7 look far apart?* Rounding. Unrounded 1.839 vs 1.738, difference 0.10, combined uncertainty 0.14 (slide 8).
- *Why does the line miss the origin?* Intercept −0.05 ± 0.03, p = 0.11, not significant. Through the origin a = 1.80 ± 0.01 (slide 8).
- *Is it linear?* An ω² term is 0.5 standard errors from zero, p = 0.65 (slide 8).
- *Why add uncertainties in quadrature?* Independent sources. Linear addition (worst case) still gives 0.6σ (slide 9).
- *Why half the range for the method?* A convention. Alternatives give 0.031 or 0.108, and the conclusion does not change (slide 9).
- *What is SEM?* Standard error of the mean, SD/√n (Part 1.5).

**Random pixels**

- *Why is your result different from the pixel method?* The rings move half a fringe, so each pixel sees less than one cycle, and fixed pixels read center and spacing changes as phase. Half the difference comes from fixed positions (slide 10).
- *The pixel mean is closer to the prediction.* The prediction cannot tell 1.70 from 1.84. The issue is the ±0.2 spread of a single draw (slide 10).
- *Would more pixels help?* They reduce the spread (0.06 at 80 pixels) but not the bias from fixed positions (slide 10).

**Drift and limitations**

- *Could the phase drift in time?* Yes, and with speeds only going up we cannot separate a steady drift from the slope. The joint fit gives D = −0.8 ± 0.5 mrad/s, p = 0.12, correlation −0.99 (slide 12, B4).
- *What causes the wander?* Plausibly the fiber's several modes. Not tested (B4).
- *What would you do differently?* Rest photos between speeds, reversed rotation, random order, count the turns, single-mode coupling with a longer-focal-length lens, record the rotation direction (slide 12).

**Speed of light**

- *Is this a measurement of c?* It is c inferred through the Sagnac formula with the given N, R, λ. The same comparison as the slope, not an independent test (slide 11).
- *What limits it?* The turn count, 7.7 % of 8.7 % (slide 11).
- *What does the slope say about N?* N = 13.8 ± 0.5 with c, R, λ known (slide 11).
- *Why is λ not in the budget?* It was given without an uncertainty. 0.16 % per nm (slide 11).

---

# Part 4. Numbers cheat sheet

Rounded as on the slides, unrounded in brackets where comparisons use them.

**Given**

| | |
|---|---|
| N | 13 ± 1 turns |
| R | 57.0 ± 0.5 cm, A = πR² = 1.02 m² |
| λ | 640 nm (no uncertainty) |
| Fiber | SMF-28e+, core radius 4.1 µm, NA 0.14, V = 5.6 at 640 nm |
| c₀ | 299 792 458 m/s (defined) |

**Prediction**

| | |
|---|---|
| $a_{\mathrm{th}} = 8\pi^2NR^2/(\lambda c)$ | 1.7 ± 0.1 rad/(rad/s) (1.738 ± 0.137, 7.9 %) |
| Δt at 1.81 rad/s | 1.07 fs, half the 2.13 fs optical period |
| Δφ at 1.81 rad/s | 3.2 rad (3.15) |

**Calibration** (slide 4, B1)

| | |
|---|---|
| Linear fit | k = 0.02039(4) rad/s per unit S, ω₀ = −0.0212(6) rad/s, χ²/dof = 47.2/7, p = 5 × 10⁻⁸ |
| Stopwatch | 0.04–0.18 % of ω, $u(T) = \mathrm{SD}/\sqrt{n(n+1)}$, lag-1 −0.34 (expected −0.43 ± 0.09) |
| Reproducibility | ε = 0.48 % of ω |
| S = 90 | timed 1.772 rad/s (8 turns, ⟨T⟩ = 3.545 s), line 1.81, quadratic 1.78 |
| Quadratic | χ²/dof = 6.0/6 (by construction), slope +1.7 % |
| Cubic | χ²/dof = 2.6/5, p = 0.75 |

**Data and rings** (slides 5–7)

| | |
|---|---|
| Photos | 161: 10 at rest, 15–22 at each of 9 speeds, 16:21:40–16:33:46 |
| First dark ring | 258 px at rest → 101 px at the top speed |
| Ring finder | 170 rays, −75° to 180°, minima within 3–8 px (rms), radii ±1 px |
| Spacing P | median 1.07 × 10⁵ px², per setting ±1.7 % |
| η | 0.62 at rest → 0.10 at the top speed (0.52 fringe, 3.3 rad) |

**Main fit** (slide 8)

| | |
|---|---|
| a (regression and calibration line) | 1.84 ± 0.02 (1.839 ± 0.024) |
| s, distance of points from the line | 0.045 rad (unrounded) |
| R² | 0.9986 |
| vs prediction | 0.7σ, p = 0.47 (z = 0.72), a/a_th = 1.06 ± 0.08 |
| Photo-error fit | χ²/dof = 56.6/8, p = 2 × 10⁻⁹ |
| Intercept | −0.05 ± 0.03 rad, p = 0.11. Through origin a = 1.80 ± 0.01 |
| ω² term | p = 0.65 |
| Wander check | 0.039(7) rad. Distance from line vs per-setting mean P: r = +0.81, p = 0.005 |

**Methods and final slope** (slide 9)

| | |
|---|---|
| Ring minima / cosine / template | 1.84 / 1.88 / 1.77 (1.839 / 1.876 / 1.768), R² 0.9986 / 0.9982 / 0.9993 |
| Method uncertainty | half range 0.054 unrounded (3.0 %) |
| Calibration model | 0.032 unrounded (1.7 %) |
| Regression and calibration line | 0.024 unrounded (1.3 %) |
| **Final a** | **1.84 ± 0.07** (0.067, 3.7 %), 0.7σ from prediction, p = 0.51 |

**Random pixels** (slide 10)

| | |
|---|---|
| 1000 draws of 10 pixels | mean 1.702(6), SD 0.196(5), fit error 0.057, ratio 3.4 |
| Middle 68 % | 1.52–1.89 |
| Failures | 2 negative (left out), 18 at the iteration limit (kept) |
| More pixels (20/40/80) | SD 0.091 / 0.073 / 0.056, fit error 0.041 / 0.040 / 0.037 |
| Own center and P (purple) | mean 1.78, SD 0.17 |
| Original single draw | 1.79 |

**Speed of light** (slide 11)

| | |
|---|---|
| c | (2.8 ± 0.2) × 10⁸ m/s (2.834 ± 0.246, 8.7 %) |
| c/c₀ | 0.95 ± 0.08, −0.7σ, p = 0.50 |
| Budget | N 7.7 %, method 3.0 %, R 1.8 %, calibration model 1.7 %, regression 1.3 % |
| Our measurement only | 3.7 % |
| N from the slope | 13.8 ± 0.5 |
| λ | 0.16 % per nm (not included) |

**Results and drift** (slide 12, B4)

| | |
|---|---|
| Largest shift | 3.29 ± 0.06 rad = 0.52 fringe at 1.81 rad/s, Δt = 1.1 fs |
| Joint fit | a = 2.1 ± 0.2, D = −0.8 ± 0.5 mrad/s, p = 0.12, ρ = −0.99, R² = 0.9991 |
| Within settings | −6.9 to +8.2 mrad/s, 6 of 10 with p < 0.05, mean −1 ± 1 |

---

# Part 5. Glossary

- **a, a_th** the measured and predicted slope of Δφ against ω, in rad/(rad/s) (= seconds).
- **Bootstrap** estimating the uncertainty of a summary number by resampling the data with replacement.
- **c₀** the defined speed of light, 299 792 458 m/s.
- **Calibration model** the choice of a straight line or a curve for ω(S). Its effect (1.7 %) is an uncertainty.
- **χ²/dof, p** the fit test against the points' own error bars (Part 1.5).
- **Fringe** one full cycle of the interference pattern, 2π of phase.
- **Gaussian beam** the beam shape of a laser and, approximately, of the light leaving our fiber: Gaussian profile, spherical wavefronts of radius R(z).
- **Gouy phase** the extra on-axis phase ζ(z) of a focused beam.
- **η** the fractional ring order at the center, φ_c = π − 2πη.
- **Method uncertainty** half the range of the three full-image slopes.
- **Mode** a field pattern a fiber can guide. Single-mode for V < 2.405.
- **P** the ring spacing in r², $2\lambda/|1/R_1 - 1/R_2|$.
- **φ_c** the phase difference at the center of the ring pattern.
- **R²** the fraction of the variation that a fitted line accounts for.
- **Regression** a least-squares straight-line fit with equal weights. Its slope error comes from the distance of the points from the line.
- **Reciprocity** both directions through the same path see the same phase. Rotation breaks it, which is the Sagnac signal.
- **rms** root mean square, here the typical distance of the minima from their fitted circle.
- **s** the regression's per-point uncertainty, 0.045 rad (unrounded), used on slide 8.
- **Residuals** the differences between the data and a fitted line or curve.
- **S** the motor controller setting.
- **SD, SEM** standard deviation and standard error of the mean, $\mathrm{SEM} = \mathrm{SD}/\sqrt n$.
- **σ (as in 0.7σ)** a difference in units of its combined standard uncertainty.
- **Splits** stopwatch lap times taken at each pass with one press, so neighboring turns share a press.
- **V number** $(2\pi/\lambda)\,a\,\mathrm{NA}$, decides how many modes a fiber guides.
- **ε** the speed reproducibility of the table, 0.48 % of ω.

---

# Part 6. Sources

- Lab manual, `Documents/Sagnac_manual.pdf`: Eq. 1–5 (Sagnac delay and phase), Eq. 6 (two-beam interference), Eq. 7 (phase from fringe shift), Eq. 9 (V number), p. 10 (optimum fiber coupling, V = 5.6 at 640 nm), p. 11 (fiber datasheet), p. 12 (setup figure, given parameters, minimum tasks).
- Preparation sheet, `Documents/Entry - Sagnac effect experiment.pdf`: Q1–Q3 (interference, Gaussian beams), Q5 (two ways to extract Δφ), Q6 (polarizers), Q7 (lens-to-fiber distance), Q8 (rotation axis), Q10 (refractive index), Q11–Q12 (fiber, single-mode at 640 nm).
- QE3 lecture slides, `Documents/QE3.pdf`: Gaussian-beam intensity (PDF p. 7) and phase (PDF p. 10).
- E. Schwartz, "Sagnac effect in an off-center rotating ring frame of reference", Eur. J. Phys. 38, 015301 (2017): Eq. 1–3 (Sagnac phase), Eq. 4 (axis independence), Eq. 5 (phase from fringe shift), Eq. 6–7 (results by linear regression), Sec. 2 (refractive index), Sec. 3 (ring-like fringes).
- G. Pascoli, C. R. Physique 18, 563 (2017): Sec. 1–2 (history, Sagnac's measurements in both directions, Fresnel drag), Sec. 3 (relativistic interpretation, Langevin, von Laue).
- A. Tartaglia and M. L. Ruggiero, Am. J. Phys. 83, 427 (2015): Eq. 1–2 (Sagnac delay, lowest order in v/c), Sec. I (gyroscopes).
- Our analysis: `RESULTS.md`, `SUMMARY.md`, `fringes.py`, `calibration.py`, `sagnac_update/`, `presentation_slides/source/` (`slides.py`, `deck.py`, `random_points_scan.py`, `collect_numbers.py`, `numbers.json`).
