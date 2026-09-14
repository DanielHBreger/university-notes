# `chua/shilnikov.py`

N6 of the plan: the eigenvalues at all three equilibria from the measured slopes, and the check |sigma| < |gamma|. Output: `shilnikov.txt` beside the script.

## What is being checked

At a fixed point the Jacobian of the three circuit equations has one real eigenvalue gamma and a complex pair sigma ± j·omega. The double scroll needs two kinds of saddle-focus: the origin with gamma > 0 and sigma < 0 (trajectories spiral in on a plane and leave along a line) and the outer equilibria with gamma < 0 and sigma > 0 (they arrive along a line and spiral out on a plane). Shilnikov's theorem says that a homoclinic orbit to a saddle-focus with |sigma|/|gamma| < 1 carries horseshoes, so the ratio is the quantity to report.

## `eigen_rows(g, Rt, rL, c1, c2, l)`

For each root of `simulate.equilibria` (the load line through the origin against the element), `simulate.stability` gives the eigenvalues and the local slope. The real eigenvalue is the one with |imaginary part| < 1 (the pair has |omega| of 10⁴ rad/s, so the split is unambiguous); the pair's real part is `sigma`, its imaginary part `omega`. `ratio` is |sigma|/|gamma|. A test checks the classification at 700 ohm and that the eigenvalues are those of `stability`.

## `kind(row)`

The saddle-focus type from the signs, in words.

## `models(identified)`

Two element/component sets. The plan's model: `simulate.make_g()` (the V–I segments) with the module's C1 = 11.5 nF, C2 = 100 nF, L = 18 mH and rL = 0. The identified circuit: `simulate.bench_static_g` with the small-signal values `bench_small_signal` returns (apparent C1, C2, L0, r0).

## `main()`

For each model and each `--r`: one row per equilibrium with v1*, the local slope, gamma, sigma, the pair's frequency in kHz, the ratio and the type; then a verdict line saying whether the geometry is the double-scroll one and whether the condition holds at the origin and at the outer equilibria. Where the load line meets the element only once (the plan's model below 403 ohm) it says so.

What the table shows: the condition holds at the outer equilibria everywhere (ratio 0.02 to 0.04), and at the origin down to about 350 ohm in the identified circuit, where gamma of the origin has shrunk enough for the ratio to pass 1. The measured double scroll ends at 328 ohm.
