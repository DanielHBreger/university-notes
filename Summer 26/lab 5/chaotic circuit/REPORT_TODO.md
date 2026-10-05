# Lab 5 report: open review items

Last reviewed: 2026-09-26, against the draft `Lab_5___Chaos (11).pdf` (sections 2.1-2.3 written).
Tick an item (`[x]`) or delete it when the draft addresses it. Numbers here match `RESULTS.md`.

Rules for the report: every number with its uncertainty and what kind it is; uncertainties
rounded to one significant figure and the value to the same digit (the user's course convention),
but every comparison computed from unrounded values; every chi-square with dof and p; figures in the report style (`report_style.py`, column width unless `figure*`).
`chua/rpot.py` is the handed-out file and is not modified, and its internals are not described.

Fixed in draft (9): 327.92(2) → 327.9(2); 2.7(5) → 2.70(5); the "352(2) Ω" threshold sentence;
the small orbit in the hysteresis paragraph; the R_pot procedure paragraph; "We measured",
"1.0029(1). with", "record-to-record deviation", "reads 0.29(1) %", "carry this factor".
Fixed in drafts (10)-(11): hysteresis range 314.4(3)-680(5) Ω; "V1 rises"; the "knob only turned
one way" sentence; the small orbit in the Fig. 3 walk-through; Table 2 doublings and chaos start;
in 2.3 Eq. (2), the period-8 claim, the measured frequency, "the a series", "discussed in Section 3".

## BLOCKING before submission

- [x] **Rewrite 2.4 in your own words** (asked 2026-09-27; done in draft 27, 28 September). You are unhappy with it and used my
      suggested text too directly; draft (14) is almost word for word the suggestion saved in
      `REPORT_TODO_2.4_suggestion.tex`. Do not submit until 2.4 has been rewritten; check the new
      draft against that file. The content can stay (M4, M5, N6, Figs. 5-6); the wording and
      structure should be yours.
- [x] **Reword 2.5 in your own words** (done in draft 30, 28 September; the C2 and inductor
      paragraphs rewritten by the user). (Added 2026-09-27, same reminder as 2.4.) Draft 27: the
      method, ideal-model, C1 and results paragraphs are the user's; the last two sentences of the
      C2 paragraph and most of the inductor paragraph are still my 28-September rewrite. Clear this
      item only when the user says so. The C2 paragraph,
      the inductor, element and results paragraphs, and the Figure 7 and Table 3 captions were
      given as text in chat (their numbers also in `REPORT_2.5_suggestion.tex`). Do not submit
      until 2.5 has been reworded; the content and numbers can stay.

## Draft (22) review (27 September; Results and Discussion merged)

- [ ] 2.5: C1's law missing, so "the growth of C1" in the results paragraph is undefined.
- [ ] 2.5: the identified model's element not stated; 2.4's "oscillator's -0.762 mS" is unexplained.
      One sentence in 2.5 (element from the oscillator records, inner slope -0.762 mS, within 1.5 %
      of the schematic) with a pointer from 2.4 fixes both.
- [ ] 2.5: Table 3 never cited in the text; Figure 7 not included (fine if dropped).
- [ ] 2.5: "The method it used" (C2) has no antecedent; C1 sentence "If we multiply ... we're left
      with" has no main clause; I_{R_t} then I_R; "As for the capacitance C1, The".
- [ ] "??" in 2.4 and Fig. 6's caption; [insert cref] / [cref to ...] placeholders in 2.1, 2.2, Table 2.
- [ ] 2.3: split done (draft 23); still "frops" and the space before ")".
- [ ] Model names: ideal / identified (as in Table 3) instead of first / 2nd / original.
- [ ] Missing space after V1, V3 ("V1is", "V1rises", ...): probably only the PDF text extraction
      (the same spot had a space in draft 22); check visually, ignore if the page looks right.

## Draft (23) review (28 September)

- [ ] 2.4: the "-0.762 mS" sentence (bound 320 Ω, Shilnikov limit ~350 Ω) was removed, so "Both of
      these limits use Table 1's G_a." now ends the discussion with no conclusion. Put it back.
- [ ] The element sentence (moved to 2.5 in draft 24, good) still says the element is "the five
      segment description from Section 2.1": wrong, it is a five-segment curve measured from the
      oscillator records, with its own slopes (inner -0.762 mS).
- [ ] 2.2 exclusion: "2 orders of magnitude smaller than the largest at 89.50(1) mV": 89.50(1) mV is
      the smallest analysed oscillation, and the 0.53(1) mV is the largest 3 kHz component of the
      excluded records ("at most").
- [ ] 2.2: "The small excess in V3 reads ..." still has no working verb.
- [ ] Figure 7: titles fixed (draft 24). Still: cite it and Table 3 in 2.5; "Local simulated
      maxima" -> "Simulated local maxima"; it floats to page 7 after the Acknowledgments (define the
      figure* block in 2.4, after Fig. 6).
- [ ] 2.1: "The op-amp in the Chua diode saturate" (draft 24) -> "The op-amps ... saturate";
      "op-amps positive" -> "op-amps'".

## Draft (43) review (3 October)

Acknowledgments now meet APS/IOP/AIP/Optica (developers, versions, what the code covered incl.
figures, how directed and verified, responsibility); onset AI note added in 2.3; 1.1 subsections
numbered (the "1.1.2" reference resolves); abstract typos, R_t, "hus", Eq. (10) brackets fixed;
Table 2 cref -> Fig. 6.
- [x] Crefs to Eq. (4) and Eq. (10) fixed in draft 44 (labels moved inside the equations).
- [ ] Table 3 still on page 9 after the Conclusions (from draft 40).
- Small (acknowledgment): no usage dates (only Springer Nature/Optica ask); "the models also carried
  out" -> as far as my records go those three analyses were Claude's; "(Which" lowercase; the 2.3
  note sits as a fragment inside the long parenthesis.

## Draft (40) review (3 October; Section 1.1 now in the report)

Abstract rewritten without the model numbers (6.6-72.1 issue gone). Intro points of 3 Oct ignored
by the user's request.
- [ ] Section 1.2 still empty; Section 2 refers to it (2.1, 2.2, 2.5) and to its two equations.
- [ ] Placeholders now resolvable: 2.1 "[insert cref]" -> Eq. (4) (eq:three_equilibria); "??" in 2.4
      and Fig. 6 -> Eq. (10) (label it eq:lyapunov_map or change the refs); Table 2 "[cref to
      lyapunov]" -> Fig. 6. 2.2's divider/rpot crefs wait for 1.2.
- [ ] 1.1 Lyapunov part: "mentioned in paragraph 1.1.2" -> the 1.1 subsections are unnumbered.
- [ ] Table 3 now floats to page 9, after the Conclusions (far from 2.5).
- Small: abstract "happend", "mismatch ... remain"; 1.1 "Rt", "INR(V1)" not subscripted, "R_t that"
  missing space, "hus"; Eq. (10) typeset with < > (use \langle \rangle, \ln).

## Section 1.1 third draft (Hila's PDF, 3 October)

Fixed: λ sign sentence; Rosenstein mentioned; stretching and folding now in 1.1.2 (the 1.1.4 reference
resolves); Shilnikov's condition |σ| < |γ| stated; 1.1.3 rewritten.
- [ ] Shilnikov part: "In the center ... it spirals outward": wrong for the origin (γ > 0, σ < 0:
      spiral in, pushed out along one direction); spiraling outward is the outer equilibria
      (σ > 0, γ < 0). 2.4 applies the condition at all three, so give both cases.
- [ ] "This orbit ... prevents the dynamics from escaping and bounds the system's motion": not what
      the homoclinic orbit does; Shilnikov: near it there are infinitely many periodic orbits
      (horseshoes), i.e. chaos. Boundedness comes from the outer saturation segments.
- Optional: one clause on what Rosenstein's method does (from the time series, separation of nearby
  points); "the diversion becomes rapid" -> the splits come faster and faster.
Still open: R_t = R0 + R_pot; "complete solution"; "two attractors" (1.1.1, end of the double-scroll
paragraph in 1.1.3); hysteresis "unstable state"; garbled inequality; equation labels.

## Section 1.1 second draft (Hila's PDF, 2 October)

Fixed: outer-saturation sentence removed; linearization mentioned. 1.1.3 being rewritten by Hila.
Still open from the first review: R_t = R0 + R_pot (and R_pot = R_Hopf - R0); "complete solution"
-> solution of the linearized equations (and linearize before the eigenvalues); "two attractors"
in 1.1.1 and 1.1.3; labels eq:three_equilibria, eq:lyapunov_map.
New in 1.1.4:
- [ ] Last sentence garbled and wrong: "If λ < 0 ... chaotic, if λ > 0" -> λ < 0 converge (periodic),
      λ > 0 diverge exponentially (chaotic).
- [ ] Hysteresis: "if we come from an unstable state" -> the coexisting attractors are the large
      cycle (periodic) and the scrolls (chaotic); a chaotic attractor is not an unstable state.
- [ ] "stretching and folding ... mentioned in paragraph 1.1.2": 1.1.2 does not mention it.
- [ ] Still missing: Rosenstein's method and Shilnikov's condition (2.4 cites both in Section 1.1).
- Optional: "pulling the trajectory back into the active inner region" (the cycle near onset stays
  on the shoulder); define λ via e^{λt} (the stray e^{λt} in the heading).

## Section 1.1 draft review (Hila's markdown, 1 October)

- [ ] Limit cycle "forced into the diode's outer saturation regions": wrong for this circuit; the
      period-1 cycles (V1 -4.7 to -0.5 V) and the double scroll never reach the outer breakpoints
      (-6.78 / +6.03 V), only the large cycle does. Growth is limited by the element's nonlinearity
      on the shoulder (its slope changes along the cycle; the curved shoulders of 2.1).
- [ ] x(t) is the solution of the equations linearized about an equilibrium, not "the complete
      solution"; linearize first, then the eigenvalues. Optional: cut the x(t) equation.
- [ ] Origin "saddle point ... between the two attractors" -> a saddle-focus, between the two scrolls
      of one attractor (ties to Shilnikov in part 4).
- [ ] Define R_t = R0 + R_pot; R_Hopf is a total resistance (R_pot = R_Hopf - R0); f = ω/2π.
- [ ] Three-equilibria inequality garbled in the markdown; label it eq:three_equilibria.
- [ ] Each breakpoint is where one op-amp saturates (inner: stage B, outer: stage A), not only Gc.
- [ ] Part 4 still to write: hysteresis (coexisting attractors, large outer cycle), return map,
      Lyapunov exponent + eq:lyapunov_map (the ?? in 2.4 and Fig. 6), Rosenstein, saddle-focus signs
      and Shilnikov |σ| < |γ|. Then 1.2 (list sent 28 Sep).

## Draft (35) review: Abstract (1 October)

Fixed in draft 37: the unit, 269.1(6) Ω, the double scroll added. Optional items left by choice.
- [ ] Draft 37: "improved the offset from measurements to 11.1(9) Ω to 58.6(5) Ω" does not match
      Table 3: the identified column (after the onset, as for the ideal range) runs from 6.6(6) Ω
      (period-3 start) to 72.1(6) Ω (large cycle, down). Use "6.6(6) Ω to 72.1(6) Ω".

## Abstract plan (28 September; content given in chat, the user writes the prose)

One paragraph, 150-200 words, no figure or section references, only numbers already in the report.
(1) what: the Chua circuit, its route to chaos as the coupling resistance changes. (2) how: I-V of
the element, then sweeps of R_pot both ways, R_pot read from the divider without opening the circuit.
(3) observed: onset 907(5) Ω, doublings to period 8, chaos from 748.3(7) Ω with a period-3 window,
double scroll 669.3(9) -> 327.9(2) Ω, hysteresis 314.4(3)-680(5) Ω; δ1 = 4(1) vs 4.669; chaos
confirmed by thin return maps and positive λ (2.4(2) to 3.1(1) x 10^3 /s). (4) models: ideal
178(2)-269.1(6) Ω high; identified (C1 10.78(3) -> 12.56(7) nF with the oscillation) onset within
17(6) Ω, doublings within 11.1(9) Ω. (5) takeaway as in the conclusions. Recheck sentence (1) once
Section 1 exists.

## Draft (33) review: Conclusions (28 September)

Fixed in draft 34: 269.1(6) Ω; "the model's onset"; δ1 rather than "Feigenbaum's constant"; "long
past" -> "past"; C1 growing with the oscillation added; the takeaway paragraph rewritten (qualitative
agreement from the ideal model, quantitative with the actual component values). Conclusions done.
Kept as is by the user (flagged once, do not raise again): the bounds sentence under the simulation
paragraph without the -0.762 mS resolution; "positive through the entire chaotic regime" (one
exception, forward 415.6 Ω, λ = -24(6) /s); the optional frequency and improvements sentences.

## Conclusions plan (28 September; content given in chat, the user writes the prose)

Two short paragraphs, ~200 words, only numbers already in Section 2. P1 what was observed: the
element (five segments, curved shoulders); the route (onset 907(5), doublings to period 8, chaos from
748.3(7) with a period-3 window, double scroll 669.3(9) -> 327.9(2), large cycle); hysteresis
314.4(3)-680(5); δ1 = 4(1); chaos confirmed (thin return maps, λ 2.4(2) to 3.1(1) x 10^3 /s,
direct method positive through the chaotic range). P2 theory and models: onset frequency within
2.8(2) % but onset 98(7) Ω off; double scroll below the 400(40) Ω bound and 420(40) Ω Shilnikov
limit (320 Ω and ~350 Ω with -0.762 mS); ideal model 178(2)-269.1(6) Ω high; identified model onset
17(6), doublings within 11.1(9), mostly from C1 growing (constant C1: up to 121.6(5)); double scroll
still ends 58.6(5) high; takeaway. Optional: improvements (fixed V/div, measure C1/C2/L first,
finer steps near the doublings). Check the draft's conclusions against this.

## Draft (31) review (28 September)

Fixed in (31): C2 now splits voltage and current; L0's origin restored. "jump" kept by choice.
- [x] Acknowledgments redone (draft 32, user's version): who did what, onset extrapolation,
      suggested text rewritten. The rpot.py-revision sentence was left out by the user's choice
      (outputs identical to 0.005 Ω); mentioned once, do not flag again.
- Optional: L0 "chosen manually" -> it was solved for numerically to match the frequency.

## Draft (30) review (28 September; physics and numbers only, as agreed)

- [ ] C2: "we split the recorded V2 into a fourier series" -- V2 alone cannot give C2; the current
      reaching node 2 has to be split too (the method compares current to voltage per frequency).
- [ ] Inductor: "grow in proportion to the jump in the current" -> the swing of the current around
      its average ("jump" reads as a sudden step).
- [ ] Inductor: L0 = 18.9(2) mH lost where it comes from; next to "rises from 18 mH" it looks
      inconsistent (zero-current value above the 0.7 mA one). Restore: chosen so the model
      oscillates at the measured frequency of the first oscillation (Section 2.3).
- Minor, user's call: "area of that circle" -> loop; low/high-frequency "matching" is a
  simplification (both are fitted to every frequency); "reducing the effects on the inductor" ->
  the inductor is closest to linear there; "fourier" -> "Fourier".

## Draft (28) review (28 September)

Fixed in (28): Fig. 6 caption reason (lobe switches up to almost twice a turn); 2.4 last sentence
(bound now fits, Shilnikov limit closer); 2.5 results softened ("improves on") with the remaining
differences. The user is ignoring the remaining Table 3 issues (do not flag), except:
- [x] Table 3 caption sign (fixed in draft 29: "Δ = model - measured").
- [x] Element sentence restored in draft 29.
- [ ] 2.5 rewording (draft 29): the user is unsure how; skeleton of facts given 28 Sep for the end
      of the C2 paragraph and the inductor paragraph, to write from with the current text closed.
- Optional, user's call: Fig. 7(b) never cited; 2.4 "for the rest of the records" (Rosenstein was run
  on every record); C2's "split ... into the harmonics" is still my wording, which the user asked
  about (simpler sentence given 28 Sep).

## Draft (27) review (28 September; 2.4 and most of 2.5 reworded by the user)

- [ ] 2.4: "which still doesn't match the measurements": wrong for the bound. With -0.762 mS the
      three-equilibria bound is 320 Ω, below the observed end 327.9(2) Ω, so it matches; only the
      Shilnikov limit (~350 Ω) is still above the end.
- [ ] 2.5 results: "fixes all three issues" while the remaining differences (double scroll end
      +58.6(5), large cycle +72.1(6) / -48(5) Ω, no small orbit) are no longer mentioned; soften
      or add one sentence. Fig. 7(b) no longer cited.
- [ ] Table 3 caption lost the sign convention (Δ = model - measured) and the unit (Ω); also the
      onsets come from linearization, not from the ±0.5 Ω steps.
- [ ] Fig. 6 caption: "(This was done because the experiment guide specifically mentions this)":
      give the reason instead, which is what the plan asks the report to state (M5 trap). Checked
      28 Sep: a turn that switches lobes takes 1.83x / 1.88x / 1.36x an ordinary turn (median) at
      650 / 549 / 450 Ω, switches are 14 / 18 / 31 % of turns, mean over median +10.5 / +16.6 /
      +1.8 %. So "up to about twice", not "about two turns" everywhere.
- [ ] 2.4: Rosenstein was computed for every forward record, not "the rest of the records".
- Kept by choice or minor (do not flag): "equilibirum", "between ... to", "slope of the loop".

## Plan M/N coverage check (draft 26, 28 September)

All of M1-M5 and N1-N6 are answered in Section 2. Gaps against the plan's wording:
- [ ] M1 intercept trap: the plan asks to check it; the report only gives -0.41(1) mA in Fig. 1's
      caption. One sentence: 0.7 of one current quantization step; a real element passes no current
      at V = 0, so it is a channel offset, removed before the curve is used.
- [ ] M4 method: how the maxima are found (Savitzky-Golay over 1/20 of a period, maxima above 2 % of
      the record's range) is stated nowhere; Figs. 3-6 all use it. One clause in 1.2 or 2.2.
- [ ] M3 record length (200 or 500 ms, 570-1400 windings): belongs in 1.2.
- Optional: M1 drive 9.17 Hz is below the plan's 20-100 Hz (harmless: the directions agree); M5's
  second trap stated without its reason (a per-branch line fit misses the turning point, where
  the slope drops toward zero); M2 "written before looking" not recorded, so do not claim it.

## Draft (26) review (28 September)

Fixed in (26): Table 2 period-3 values (correct as start/end in each sweep's own direction: fwd
718.9 then 713.1, back 707.1 then 712); Table 3 period-3 rows; Figure 7 placed before 2.5 and its
caption; element sentence now gives -0.762 mS; Acknowledgments "with them". (My note on draft 25
that 718.9 and 712 must share a row was wrong for start/end rows.)
Still open:
- [ ] "??" twice and the cref placeholders: resolve once 1.1 and 1.2 are written (user, 28 Sep).
Kept as is by the user's choice (28 September; minor, or reads human): do not flag again.
- Table 2 labels "Period 3 up / down"; -0.762 mS without an uncertainty; the 2.5 element sentence's
  wording; 2.1 "The op-amp ... saturate", "op-amps positive"; 2.2 "The small excess in V3 reads";
  2.3 "frops", "Section 2.5 )"; 2.5 "(Table 3)..", "at a 178(2) Ω", "The method it used";
  Acknowledgments "it also carried out", "It reviewed".

## Draft (25) review (28 September)

Fixed in (25): 2.4's -0.762 mS sentence restored; C1 half-sentence in 2.5; Fig. 7 and Table 3 cited;
2.2 exclusion now "smallest".
- [ ] Table 2: rows renamed "Period 3 up / down" and the FORWARD column swapped too (now 713.1 in "up",
      718.9 in "down"). Correct: upper edge 718.9(3) (fwd) and 712(2) (back); lower edge 713.1(5)
      and 707.1(5). Label "upper edge" / "lower edge" ("up/down" reads as sweep direction).
- [ ] 2.4 now points to Section 2.5 for -0.762 mS, but 2.5 never gives it: fix the element sentence.
- [ ] 2.1: "The op-amp in the Chua diode saturate, so the shoulders curve" lost "softly", which is
      the explanation: "The op-amps in the Chua diode appear to saturate gradually, so ...".
- [ ] 2.5: "(Table 3).." double period.

## Draft (24) review (28 September)

- [ ] Table 3: the period-3 rows were swapped in draft 24 (now "start 713.1(5)", "end 718.9(3)"),
      which contradicts Table 2's forward column. Restore 718.9 first; label both tables' rows
      "upper edge" / "lower edge" so the direction cannot confuse. The back column of Table 2 is
      still the one that is swapped (going up: 707.1(5) first, then 712(2)).

## AI disclosure (asked 2026-09-27)

Acknowledgments text given in chat. In-text mentions suggested (check each draft):
- [x] 2.3 period-8 check: "An AI assisted analysis of the records (see Acknowledgments)" (draft 40).
- [x] 2.3 the onset from the A1^2 extrapolation and its ±5 Ω calibration systematic (onset.py, Claude)
      (in draft 52, inside the chi2 parenthesis).
- [x] 2.5 "The component values were identified from our records with an AI assistant" (draft 40).
- Checked against journal policies on 3 Oct (APS June-2026 update, AIP, IOP, Springer Nature,
  Elsevier, Science, Optica): add developers and usage dates, list what the code covered (incl. the
  figures), state how outputs were directed and verified (APS, IOP, Optica); research-process uses
  where the method is described (onset in 2.3; optional: the simulation code in 2.5). Text in chat.
- [ ] Everything else (2.1 fit code, error bars, figures, Table 2's transitions, 2.4's analyses, the
      simulation, draft review and suggested text) covered by the Acknowledgments. The git history
      also shows code from a ChatGPT review (9 September): the user decides whether to name it.
- [x] Resolved in draft 54: Acknowledgments now "edited or rewrote in our words"; the "like the
      swing" reason and "as the loop area does" restored.
- 5 Oct (inductor law, Eq. (23) in draft 52): no separate disclaimer needed. 2.5's sentence is
      now "The effective parameters were identified ... with an AI assistant", which covers the
      laws too. Open: the paragraph after Eq. (23) is my rewording verbatim, but the
      Acknowledgments say suggested wording was "rewrote in our own words". Either reword it or
      amend that sentence (text given in chat 5 Oct).
      Draft 53: light word-level edits only (same for the Eq. (22) lead-in); Acknowledgments
      unchanged, so still mismatched. Suggested "which we reviewed and then edited or rewrote".
      Draft 53 also dropped "like the swing", so "so d is measured from the running average"
      no longer follows; "grows with the current swing, with the loop's area" reads oddly.

## Draft (45) review (3 October; Section 1 complete, Table 4 now placed in 2.5)

Section 1 added Fig. 1 and Table 1, so every later figure and table number moved up by one.
- [x] 2.2 crefs: fixed in draft 47 (R_pot = R0 B/A numbered Eq. (13) in 1.2; 2.2 cites Eq. (12)
      and Eq. (13)). No critical issues left in Section 2 onward as of draft 47.
- [x] 2.4: "Rosenstein's method (Section 1.1.3)" -> Section 1.1.4 (fixed in draft 46).
- [x] Fig. 8(a) panel title said "element of Table 1" (now the component table): relabelled
      "Table 2" in report_figures.py and simulation.pdf regenerated (uploaded in draft 46).
- [x] 1.2 said the I-V drive was 20-100 Hz (draft 46: frequency removed from 1.2; 2.1 keeps 9.17 Hz).
- [x] Conclusions overclaim (draft 46: "eliminated most of the difference from the measurements").
- Uncertainty review items 3-5 (Fig. 7 "statistical", C1 SE vs SD, 3045.06(3) Hz) not applied;
  the user's call, and only needed if the appendix is written.

## Draft (50) review (3 October): judged submittable

References now Kennedy 1993 I/II, guide, Strogatz, plan, Rosenstein, Feigenbaum; Fig. 1 "taken
from [5]"; "10 % to 20 %" cited to the plan; Shilnikov cited to Kennedy II. Nothing wrong in the
results. Recommended before submitting: crop Fig. 1's caption strip; RK4 step, Rosenstein
parameters, Table 4 (Ω); 2.4 "experiment guide said" -> "experiment plan [5]". Rest optional.

## Draft (48) review (3 October): references added (6 entries)

- [ ] Fig. 1: still no source ("Adapted from [2]"); image still carries the guide's caption strip
      (with CH1/CH2 swapped against the data) and R0 = 990 Ω vs Table 1's 992 Ω.
- [ ] Shilnikov's theorem (1.1.3, 2.4) cited nowhere: add shilnikov / silva.
- [ ] "10 % to 20 % lower" is from the plan: cite [4] in 1.1.4; 2.4 "the experiment guide said"
      -> "the experiment plan [4]".
- [ ] Not yet done from the ChatGPT list: RK4 step, Rosenstein parameters, Table 4 (Ω) headers,
      the 1.1 wording fixes, next-step sentence.
- Cosmetic: "circuit[1][2]" -> "circuit~\cite{matsumoto,guide}" ([1, 2] with a space);
  "4.669[6]" needs a space; Strogatz prints "2nd" without "ed." (edition = {2nd ed.});
  ref [4] is labelled "Lab guide" (should say experiment plan).

## ChatGPT review of draft 47 (3 October): what holds

- [ ] References section + citations (none yet). Bibliography given in chat: guide, plan, Matsumoto
      1984, Chua-Komuro-Matsumoto 1986, Kennedy 1992, Feigenbaum 1978, Shil'nikov 1965, Silva 1993,
      Rosenstein 1993, Strogatz 2015, Savitzky-Golay 1964. Fig. 1 needs "adapted from [guide]"; the
      image still carries the guide's own caption ("Fig. 2. Chua oscillator...", "N_R Fig. 1"):
      crop it. 2.4's "10 % to 20 %" is from the plan, not the guide.
- [ ] 2.5: RK4 step 0.5 us (ideal), 0.1 us (identified); checked against DOP853 (rtol 1e-10) and
      a 4x smaller step over 1 ms: max |dV1| <= 0.74 mV (verify_models.txt).
- [ ] 2.4: Rosenstein parameters: 3-D embedding (V1, V2, V1(t - T/4)), unit variance, SG T/20,
      Theiler window T, pairs followed 3.5T, slope fitted 0.5T-2.5T, ~50 000 references,
      u = SEM over 4 time blocks; T = median spacing of the maxima of the record.
- [ ] Table 4: (Ω) on both Δ headers.
- [ ] 1.1 (Hila): "complete solution" -> local solution of the linearized system; "averaging these
      local slopes" -> averaging ln|f'|; 1.1.4 "unstable state" -> chaotic state.
- [ ] Optional: hidelinks; one "next step" sentence (other days' sweeps as an independent test).
- Wrong: "Table 2 cites Fig. 2, currently Fig. 4": Fig. 2 is the I-V plot in draft 47.

## Uncertainty review against draft 44 (3 October; for the planned appendix)

Checked every quoted uncertainty against the code and re-ran the ones the code does not store
(scratch scripts; nothing in the pipeline changed). Most are right. To fix:
- [x] (draft 45) Double scroll end 327.9(2) -> 327.9(3): its bracket (trace402 327.83, trace403 327.91) is
      only 0.04 Ω, smaller than the two records' own u (0.40, 0.33 Ω; they even read in inverted
      order). Include the endpoints: u^2 = a^2 + (u1^2 + u2^2)/4. Same rule: large cycle
      314.4(3) -> 314.4(4). Table 3: 255.6(5) -> 255.6(6), 58.6(5) -> 58.6(6); 72.1(6) unchanged.
      The other brackets don't change.
- [x] (draft 45) Shilnikov at the origin used the rounded G_a = -0.72 mS. Unrounded: 0.349 (+0.038/-0.031)
      -> 0.35(4); the ratio reaches 1 at 427 (+41/-39) Ω -> 430(40) Ω (2.4 and Conclusions).
- [ ] Fig. 6 caption "Error bars are statistical": the return-map bars also hold the fit-window
      term, which dominates at 549 and 450 Ω (151 of 175 and 415 of 429 s^-1).
- [ ] 2.5: 10.78(3) nF is the SE over the 4 onset cycles; 12.56(7) nF is the SD over the 48
      large-cycle records. Say which, or quote the SE for both. Both leave out the ~1 % range-to-range
      floor (~0.1 nF).
- [ ] 3045.06(3) Hz is statistical only; the timebase (25 ppm, calibration age unknown) is at
      least 0.04 Hz. Say "statistical", or write 3045.06(5). The 2.8(2) % doesn't change.
- [ ] Appendix: state that bracketed values are midpoint ± half-width (a bound; the standard
      uncertainty is half-width/sqrt(3)). δ1 holds either way: 0.9 half-widths or 1.6 standard
      uncertainties from 4.669; Monte Carlo 95 % interval [2.91, 5.15], P(δ1 >= 4.669) = 0.10.
- [ ] Optional: the left inner breakpoint -1.09(4) V is really ±0.13 V (its χ² profile is flat over
      two codes); the other three are within one code, so ±0.04 V holds for them.

## Errors (wrong as written)

- [x] **2.3 predicted onset** (found 2026-09-27; fixed in draft 23): "R_pot = 1005(5) Ω with a frequency of 2961(7) Hz"
      was computed with the rounded G_b = -0.418 mS. From the unrounded -0.41756 mS it is
      **1006(5) Ω and 2963(7) Hz**. The differences quoted after it, 98(7) Ω and 2.8(2) %, are
      unchanged.

- [x] **2.2 hysteresis paragraph** (new in draft 10; fixed in draft 23): "Turning the resistance down, the double scroll
      remains to 314.4(3) Ω, followed by a small orbit" → the double scroll remains to 327.9(2) Ω; the
      small orbit lasts from there to 314.4(3) Ω, below which the large cycle takes over. Optional:
      "a loop 366(5) Ω wide, plus the 0.5–2 % range-to-range scale (the two ends are on different
      ranges)".
- [ ] **2.2 exclusion evidence** (reworded in draft 10, still wrong): "the oscillations of V1 are 2
      orders of magnitude smaller than the value of V1". Against V1 itself (−3.3 to −4.2 V) the 3 kHz
      component (at most 0.53(1) mV) is nearly 4 orders smaller, and V1's value is not the right
      comparison anyway. Compare with a real oscillation: "the 3 kHz component of V1 is at most
      0.53(1) mV, against 89.50(1) mV in the smallest analysed oscillation" (a factor of about 170).
- [x] **2.2 scale shift** (fixed in draft 23): "R_pot jumps upward at each range change" → it jumps **in either
      direction** (+3.8 Ω at records 36/37, +9.3 Ω at 42/43, −17.6 Ω at 43/44). "its V1 indeed places
      it" → "its V1 maximum places it between records 42 and 44 (852.9(2) and 844.5(2) Ω)".
- [ ] **Table 2**: "Period 3 start / end" is reversed for the back column (turning up, the window
      starts at 707.1(5) and ends at 712(2) Ω). Label the rows "Period-3 window, upper edge" and
      "lower edge".
- [x] **2.2 first paragraph** (fixed in draft 23): "fewer measurements in the increasing direction as the qualitative
      behavior is different due to hysteresis" is not the reason. Nothing changed along the large
      cycle (3–680 Ω), so it was recorded sparsely there, and densely only between 685 and 780 Ω.
      Move to the hysteresis paragraph.

## 2.2: structure (agreed order, status after draft 9)

- [ ] (1) sweeps: done, except moving the fewer-back-records sentence to (6). Optionally add
      "apart from one short reversal near 775 Ω (Section 2.3)" after "turning the knob [...]".
- [x] (2) procedure and exclusion: written (fix the evidence sentence above).
- [ ] (2) still missing: "for the remaining 479 records the residual is below 1 % in 460 cases, and
      the statistical uncertainty of R_pot has a median of 0.2 Ω (range 0.04–2.5 Ω)". Optional:
      "R0 = 992 Ω (measured on the bench; its uncertainty was not recorded)".
- [ ] (3) A + B paragraph:
      - "Over these records" follows the sentence about the 3 weakly oscillating records → "Over
        the 479 analysed records".
      - "The small excess the channel measuring V3 reads [...]" has no verb → "The small excess is
        what a V3 channel reading 0.29(1) % higher than the other two would produce; since A and
        B both carry that factor, it cancels in R_pot."
      - Add the limitation: the check cannot detect a gain difference between the V1 and V2
        channels, which is what sets the resistance scale (leads into (4)).
      - Add the near-onset sentence: "On the 12 smallest analysed oscillations (records 14–23
        forward, 73–74 back) the resistance obtained this way reads low, by up to 13 Ω; for these
        records it is taken from the mean voltages instead (Section 2.3)."
- [ ] (4) scale shift: still the tail of the Fig. 3 paragraph. Make it its own paragraph after (3),
      fix the errors above, and add: the size (0.5–2 % between ranges; values on the same range are
      unaffected) and the forward/back comparison at R1 (different ranges, ADC steps 26.5 and
      24.1 mV; 776(2) against 768.1(2) Ω, a difference of 8(2) Ω).
- [ ] (6) hysteresis paragraph:
      - add that the back sweep never re-enters the double scroll, so the double scroll's upper edge
        on the way up could not be measured;
      - "stays on the large circle far beyond that point" → "stays on the large cycle past
        314.4(3) Ω, where it had appeared on the way down, up to 680(5) Ω";
      - "The system decides which of the two attractors to follow depending on [...]" → "Which of
        the two attractors the circuit follows depends on the direction from which R_pot is
        approached".
- [ ] Paste the `eq:divider` / `eq:rpot` block into 1.2 so the crefs resolve.

## Table 2 caption

- [ ] Say what the parentheses are: half the gap between the two records that bracket each
      transition; the onset is extrapolated (Section 2.3) and includes the ±5 Ω calibration
      systematic; comparisons across oscilloscope ranges carry an extra 0.5–2 %.
- [ ] Say what "–" means (not resolved or not sampled). Back column of "Double scroll start" →
      "not reached".
- [ ] Row "Small orbit start" → "Double scroll → small orbit".

## Figure captions

- [ ] Fig. 3: "V1in" is missing a space. "The scatter pattern near 860 Ω comes from the shift in
      the oscilloscope scale shifting between different ranges" says "shift" twice → "The scatter
      near 860 Ω comes from shifts of the resistance scale between oscilloscope ranges
      (Section 2.2), not from the circuit." Mention that the records nearest the onset are placed
      at their mean-voltage resistances (Section 2.3).
- [ ] Fig. 2: "40ms" → "40 ms".

## 2.2: typos and wording

- [ ] "Fig. 2 Shows" → "shows"; ". at 669.3(9) Ω" → capital A; "then expanding" → "and then expands".
- [ ] "maximas" → "maxima" (already plural); "a large cycle behavior" → "a large cycle".
- [ ] "an island of stability producing a narrow window of a period 3 cycle" → "a narrow period-3
      window".
- [ ] "large circle" → "large cycle"; "do not trace over each other perfectly" → "do not retrace
      each other".
- [ ] "V1falls" is missing a space (the sentence changes anyway, see Errors).
- [ ] Fill in the cref placeholders ([cref to the divider equation], [cref to the rpot equation],
      [cref to lyapunov], [insert cref] in 2.1).

## 2.1: still open

- [ ] Fig. 1 caption: "the estimated voltage quantization" → "the voltage quantization projected
      through the local slope"; add "Systematic errors (oscilloscope gain, shunt tolerance) are
      not included."
- [ ] Fig. 1 caption: move the −0.41(1) mA intercept sentence into the text with its meaning: a
      real element passes no current at V = 0, so it is a channel offset (0.7 of one ADC level),
      removed before use.
- [ ] "rising to −0.47 mS" → "steepening to −0.47 mS".
- [ ] Breakpoints in the text lack units: "(−1.09(4) vs 0.58(4))" → "(−1.09(4) V and 0.58(4) V)";
      same for the outer pair.
- [ ] Table 1 caption: "…at each end of the trace, where the drive turns around"; add "Slope
      uncertainties are statistical; on the shoulders the slope varies along the segment (see text)."
- [ ] M2, after Eq. (1): the upper bound lies beyond the 1 kΩ dial, so only the lower bound is
      testable; its uncertainty comes from G_a. Fill "[insert cref]".
- [ ] Upright subscript `R_{\mathrm{pot}}` everywhere (Eq. 1, captions, Table 2 caption). It is
      still italic.

## Moved to Discussion (check they appear there)

- [ ] M2 against the observation: double scroll 669.3(9) down to 327.9(2) Ω; the lower end is
      72(36) Ω below the predicted bound 400(40) Ω (2σ; 65(36) Ω after removing the +2 % step below
      467 Ω). Resolution: with the oscillator's own G_a = −0.762 mS the bound is 320 Ω. Same for N6 (2.4):
      with the measured G_a, Shilnikov's condition at the origin holds only down to 420(40) Ω.
- [ ] Cause of the scale shift: channel gains differ between vertical ranges (scope spec ±2 % of
      full scale).
- [ ] Why two attractors coexist (the hysteresis).
- [ ] Shoulder curvature (soft op-amp saturation) and the asymmetric breakpoints (unequal op-amp
      rails).
- [ ] G_a is the least-determined slope (13 points, ±3 %): trace −0.72(2) mS against −0.762 mS from
      the oscillator and −0.758 mS from the schematic.
- [ ] Onset against the models: measured 907(5) Ω; plan's model 1006(5) Ω; identified circuit
      925(3) Ω (2.5). (The earlier 905 Ω and the C1 = 11.5 nF model's 887 Ω are superseded.)
- [ ] From 2.3: see "Discussion items from 2.3" in the 2.3 plan below.

## 2.3: review of draft (12) (2026-09-26)

Fixed in (12): 3045.06(3) Hz; 0.9 uncertainties; the sentence ending; Fig. 4 is the current file;
"In two cases", "phase of the splitting", "splits of the V1 maxima".
Still open:
- [ ] "the squared amplitude of V1 ... reaching zero at 907(5) Ω (...), and has a frequency of
      3045.06(3) Hz": the amplitude is the subject; the first oscillation has the frequency.
- [ ] Why the onset is extrapolated (a resting circuit gives no R_pot), that the fit uses the four
      smallest oscillations, and what the ±5 Ω systematic is (calibrating R_pot for these records).
- [ ] After "δ = 4.669.": consistent, but one ratio known to ~25 % cannot test universality (the
      constant is the limit of the ratios); where Table 2's "chaos begins 748.3(7) Ω" comes from
      (between the last steady period-8 record and the first chaotic one). Beyond the guide, so
      optional under the guide-only rule (2026-09-26): the accumulation point 749(1) Ω and the
      pointer to a Lyapunov exponent at 746.2 Ω.
- [ ] Cite the Hopf formula and the definition of δ (1.1); "component values" → "nominal component
      values".
- [ ] "the last of the previous range and the first of the new range" → "the last record of the old
      period and the first of the new" ("range" reads as the oscilloscope range); "a series of period
      doublings for the oscillations' period occur" → "the oscillation undergoes a series of period
      doublings".
- [ ] Fig. 4 caption: "panel b" → "Panel (b)"; "half the summed gap ..., which introduces
      uncertainty as for the exact value" → "The shaded bands span the gap between the two records
      that bracket each transition"; "capture range" → "vertical range"; add "Each column of points
      is one record."

## 2.3 Route to chaos (consolidated plan, 2026-09-25; full explanation given in chat)

Format: three paragraphs (P2 may split in two), one column-width figure (the cascade, which becomes
Fig. 4 if the onset figure is dropped), no new table (Table 2 has the transition values; a small
predicted/measured table for N1 is optional). Numbers: half-brackets for R_n (as in Table 2),
statistical (and systematic where stated) for the rest.

Opening sentence: follow the circuit from rest to chaos as R_pot is lowered: the oscillation is
born (Hopf), doubles its period three times, and becomes chaotic.

P1 Onset of oscillation (N1). Purpose: test the linear stability prediction.
- Prediction (1.1 formula; Table 1's G_b of the left shoulder, since the circuit starts
  oscillating about the negative equilibrium; nominal C1 = 10 nF, C2 = 100 nF, L = 18 mH;
  R0 = 992 Ω): R_pot = 1005(5) Ω, f = 2961(7) Hz ((·) from G_b only).
- Measured: resting records have no R_pot, so the onset is approached from the oscillating side:
  the first oscillation is at 904.9 Ω, and its amplitude extrapolates to zero at 907(5) Ω (the
  √ law of a supercritical Hopf; 5 Ω calibration systematic). Frequency of that first oscillation
  3045.06(3) Hz (statistical; changes < 2 Hz extrapolated to the onset).
- Comparison: frequency 2.8 % above, onset 98 Ω (10 %) below.
- One sentence: the onset, unlike the frequency, depends sensitively on the element's slope at the
  equilibrium, which varies along the shoulder (2.1): each 0.01 mS moves the prediction 25 Ω and
  the frequency only 30 Hz. Details → Discussion.

P2 Period doubling (N5). Purpose: locate R1–R3 and show the cascade is real and reversible.
- Fig. 4: one maximum → 2 → 4 → 8 → bands.
- Period of a record = cycles after which the V1 maxima repeat to within the noise. A record counts
  only if the pattern is steady through the record; records 67 and 94 flip the phase of their
  splitting (too close to a doubling to place) and are left out.
- Each R_n = midpoint between the last steady record of the old period and the first of the new,
  ± half the gap (`chua/delta_bench.py`): R1 = 775.7(19), R2 = 755.6(8), R3 = 750.3(8) Ω.
- Records 63–66: knob back up to 783.0 Ω → period 1 again → steady period 2 again at 773.58 Ω
  (record 68): no hysteresis at the doubling.
- Records 61–100 on one oscilloscope range: the scale shift does not enter the differences.

P3 Feigenbaum ratio and the start of chaos. Purpose: test the universal scaling, locate chaos.
- δ1 = (R1 − R2)/(R2 − R3) = 3.8(10) (half-brackets propagated, R2 shared); the universal 4.669
  is 0.9 of the uncertainty away. δ1 is the first ratio only. (Earlier: 3.7(8) and the
  presentation's 4 ± 1, both of which used records 67 and 94.)
- Accumulation point with δ = 4.669: R3 − (R2 − R3)/(δ − 1) = 748.9(11) Ω; chaos starts at
  748.3(7) Ω (between record 93, the last steady period 8, and record 97, the first chaotic
  record), 0.5(12) Ω below. Period 16 not resolved.
- Period 8 is held steadily by three records only (92, 93, 95). Record 94 has the pattern but it
  collapses and re-forms with the opposite phase; record 96 is not steady. Evidence:
  `chua/period8_check.png` and `.txt` (χ²/dof = 725–1186 against period 4, p < 1e-300; nothing at
  16, p = 0.59–0.99; the same splittings in both halves of each record).
- Optional: bands merge 4 → 2 below that (records 97–98, 99 on), visible in Fig. 4(b).
- Back sweep: same sequence, R1 = 768.1(2) Ω on its own range (2.2).
- Pointer: the return-map Lyapunov exponent is positive from 746.2 Ω, 549(86) s⁻¹ (2.4).

Fig. 4 caption (suggested): "Local maxima of V₁ through the period-doubling cascade, forward
sweep: (a) overview, (b) enlargement. Each column of points is one record. Dashed lines mark R₁,
R₂, R₃ and the start of chaos (Table 2); shaded bands show half the gap between the bracketing
records.
All records shown were taken on one oscilloscope range."

1.1 must supply: linearisation about an outer equilibrium, with characteristic polynomial
λ³ + a2 λ² + a1 λ + a0 and, for R = R0 + R_pot and G = 1/R,
a2 = (G + Gb)/C1 + G/C2, a1 = G Gb/(C1 C2) + 1/(L C2), a0 = (G + Gb)/(C1 C2 L). The Hopf point is
a2 a1 = a0, which gives the plan's R_Hopf, with frequency ω² = a1 there. The √ amplitude law of a
supercritical Hopf (one line); period doubling, the √ growth of the splitting, δ (definition,
4.669) and the geometric accumulation R∞ = R_n − (R_{n−1} − R_n)/(δ − 1).

Discussion items from 2.3:
- Onset sensitivity. The slope at the equilibrium (V1 = −3.30 V) is uncertain by ±10–15 %:
  −0.390(10) mS (middle third), −0.33 to −0.35 (local lines), −0.429(4) (quadratic). That spans
  977–1074 Ω for the prediction or wider, with f 2924–3048 Hz. C1 moves it 8 Ω per % (C1 = 10.9–11.8
  nF would reproduce the onset, depending on the slope). The frequency is set mainly by L and C2 and
  is the robust check; the onset tests the inputs more than the theory.
- Optional: the amplitude grows gradually over about 60 Ω, where a piecewise-linear element starts
  abruptly at finite size, which points to smooth curvature (2.1).
- δ1 is the first ratio only; the models give no δ (period 4 goes straight to chaos).

## 2.4 Chaotic regime (plan revised 2026-09-26: only what the guide asks, M4, M5, N6; full text given in chat)

Figures (in `report_figures/`, made by `chua/report_figures.py`):
- Fig. 5 `return_maps.pdf` (figure*, 2 x 4): (a)-(e) one record per regime (the user's addition);
  (f)-(h) the three long double-scroll records 650.3(2), 549.3(3), 450.2(3) Ω (trace229, 270, 342),
  which carry M4 and M5.
- Fig. 6 `lyapunov.pdf` (column), the user's choice (like the presentation's Lyapunov figure):
  the direct (Rosenstein) exponent of every forward record, plus the three return-map estimates
  as open squares. The caption carries the two M5 pitfalls.
P1 M4: return maps; periodic records give points, chaotic ones narrow curves: a one-dimensional map
  whose stretching and folding make the chaos (theory in 1.1). Over (f)-(h): positive-lobe maxima
  1.84-2.50 → 1.15-1.51 V; lobe changes in 14(1), 18(1), 31(1) % of turns.
P2 M5 (two sentences + Fig. 6): direct exponent zero on the large cycle, within about
  0.3 × 10³ s⁻¹ of zero on the periodic records, rising from zero where chaos begins to about
  1.9 × 10³ s⁻¹ at the start of the double scroll, about 1-2 × 10³ s⁻¹ through it (a few dips), zero
  where it ends. Map estimates 3.1(1), 2.4(2), 2.4(4) × 10³ s⁻¹ at 650, 549, 450 Ω, 1.2(2)-1.72(7)
  times the direct values, more than the 10-20 % usually found (cite the plan) → Discussion.
P3 N6: saddle-focus signs at the three equilibria; |σ|/|γ| = 0.045-0.054 at the outer ones; at the
  origin 0.34(4) at the upper edge of the double scroll, reaching 1 at 420(40) Ω; the double scroll
  continues to 327.9(2) Ω, past this and past the 400(40) Ω three-equilibria bound → Discussion.
Wording: the user asked for less of the guide's own phrasing; the suggested text paraphrases it.
Rewritten 2026-09-26 in the report's own LaTeX conventions (main.tex: \cref, \qty and \num with \pm,
$R_{pot}$, labels VAC/sweeps/route/theory/discussion, tab:nonlinear_vac_data, tab:regime walks).
- [ ] Needs \label{chaotic} on the subsection and \label{eq:lyapunov_map} in 1.1; Table 2 caption
      "[cref to lyapunov]" → \cref{fig:lyapunov}.
Removed as beyond the guide (analyses stay in RESULTS): λ⟨T⟩, the period-3 window's 3 → 6
  doubling, the Table 2 "Period 3 → 6" row. (The exponent across the sweep is back, at the user's
  request.)
- [ ] Discussion: map estimates 1.2(2)-1.72(7) times the direct values, beyond the 10-20 % the plan
      expects (a simulated test with this circuit's parameters gives 1.14(11)).
- [ ] 1.1 needs: the return map (points on a curve = a one-dimensional map; stretching and folding,
      the horseshoe); the Lyapunov exponent and λ from the map, λ = ⟨ln|f'(M_n)|⟩ / ⟨T⟩ (label
      eq:lyapunov_map); the equilibria of the piecewise-linear circuit, their eigenvalues γ and
      σ ± iω, saddle-foci, and Shilnikov's condition |σ| < |γ|.
- [ ] 2.2 needs: how the maxima of V1 are found (Savitzky–Golay over 1/20 of a period, maxima above
      2 % of the record's range), since Figs. 3–5 all use them.
- [ ] Table 2 caption: "with chaos confirmed by the Lyapunov exponent (see [cref to lyapunov])":
      supported by Fig. 6 (the exponent across the sweep, back in at the user's request); point the
      cref to Fig. 6.

## 2.5 Numerical simulation (corrected B, 2026-09-27)

Not in the experiment plan; the user's own addition. The user chose B (the plan's model plus the
circuit with its measured parts), asked for everything under it to be verified, then to use the
corrected B. Verification record: RESULTS.md, "Simulation: verification of 27 September";
`chua/verify_models.py`. B as first described rested on errors (inductor law charging the DC
current with core loss; C1 averaged over amplitudes; L0 from a biased resistance; two op-amp
constants set against the transitions); its old agreement (905, 776.8, 697 → 344, 683 Ω) is void,
including where the presentation quotes it.

Suggested content: given in chat on 27 September (evening), revised for node 1's law;
`REPORT_2.5_suggestion.tex` is the earlier version (constant C1) and its numbers for model (b)
are superseded. Figure: `report_figures/simulation.pdf` (figure*; define it in 2.4 after the
Lyapunov figure so that it lands on the next page; regenerated with node 1's law). Numbers (Ω):

Model (b) in the report: node 1's capacitor law with the op-amps at their datasheet speed (the
user chose on 27 September to leave stage B's measured lag out of the paper; that version stays in
RESULTS.md). Table 3 as deviations from the measurement (structure agreed 27 September: P2 and P3
describe the models, Table 3 holds every number, P4 only interprets). Model minus measured (Ω):

| | measured | (a) | (b) |
|---|---|---|---|
| onset | 907(5) | +98(7) | +17(6) |
| period 1→2 | 776(2) | +178(2) | -10(2) |
| period 2→4 | 755.6(8) | +182.9(9) | -11.1(9) |
| period 4→8 | 750.3(8) | +184.2(9) | -10.8(9) |
| chaos begins | 748.3(7) | +185.2(9) | -9.8(9) |
| period-3 window, upper / lower edge | 718.9(3) / 713.1(5) | none | +6.6(6) / +10.4(7) |
| double scroll start | 669.3(9) | +217(1) | +25(1) |
| double scroll end | 327.9(2) | +255.6(5) | +58.6(5) |
| large cycle begins (down) | 314.4(3) | +269.1(6) | +72.1(6) |
| large cycle ends (up) | 680(5) | +248(5) | -48(5) |

Node 1's law (RESULTS.md, "Node 1's capacitance law"): C1 = 10.78(3) nF on the onset cycles (the
earlier 10.80(2) nF was an estimator bias), growing with the swing since the last turning point
(0.843(9) nF/V to 0.5 V, 0.1775(12) nF/V beyond: Rayleigh law of a ceramic dielectric).

- [ ] 2.5 written in the user's own words (the 2.4 lesson: do not paste the suggestion verbatim).
- [ ] Draft (18) review of 2.5 (27 September): method, model (a), C1 and C2 paragraphs in. Missing: C1's
      law (the model uses it, not the per-record values); the chaotic-record sentence (no closed cycle,
      the element's term negligible over a long record); the inductor paragraph (dropped in 18; text
      given in chat); the element sentence; the model (b) results paragraph (text given in chat); the
      figure* block (define in 2.4 after Fig. 6) and its crefs; Table 3's caption (text given in chat)
      and units in its Delta headers; one name per model (ideal/identified vs (a)/(b) vs first/2nd).
      Fixes: "I_{R_t}" then "I_R"; the "If we multiply ... we're left with" sentence has no main
      clause; 11.6(2) and 12.56(7) nF are mean and SD over records (say so); "happens at a 178(2)";
      "at 1 Ω steps" for the missing period-3 window; "remembering the direction ..." → the state is
      carried from each step to the next.
- [ ] 1.1 needs the circuit equations (the 2.5 text has "[cref to the circuit equations]").
- [ ] Discussion: with node 1's law the gaps of model (b) have both signs (-11 to +59 Ω, -48(5) at
      the large cycle's end going up) and no longer grow with the amplitude (with C1 constant:
      up to +122 Ω). Model (b)'s node 1 matches the records on the period-1 cycles but falls
      0.1-0.4 nF (1-3.5 %) short where the orbit spends time inside the inner breakpoints
      (double scroll, large cycle): the records show extra capacitance there that the model
      leaves out (it is the op-amp's delay, kept out of the paper by choice), likely part of why
      the double scroll ends 59 Ω late. C1's type is not recorded: check the part (a class-2
      ceramic would explain the law). Model (a) depends strongly on the inner breakpoints, the
      least certain entries of Table 1.

## Coming next

- Abstract, 1.1 and 1.2 are empty.
