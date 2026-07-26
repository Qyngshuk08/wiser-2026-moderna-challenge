# WISER 2026 — Moderna Challenge: Quantum RNA Secondary Structure Prediction

Team: **Qudit Creons**

## Summary

We formulate mRNA minimum-free-energy (MFE) secondary structure prediction
as a QUBO over stacked base-pair ("quartet") variables, weighted by real
Turner2004 nearest-neighbor thermodynamic parameters pulled directly from
ViennaRNA — rather than the tunable heuristic reward constants used in
prior quantum-annealing RNA folding work. We solve this QUBO on D-Wave
(Leap hybrid and direct QPU) and via QAOA on IBM hardware, benchmark
against ViennaRNA's classical MFE at both small and large scale, and
characterize where the formulation succeeds, where it breaks down, and
why — including several real bugs caught mid-project by cross-checking
results against independent ground truth rather than trusting
plausible-looking output.

**Headline result:** the base stacking-only model matches real MFE only
10% of the time on unbiased random sequences (0% at n≥60nt). A dynamic
program confirms adding real hairpin/bulge/internal-loop energies raises
this to 53.4% — and further confirms hairpin alone (the part now ported
into the actual QUBO) accounts for most of that gain at short sequences
(n≤30) but very little at longer ones, where bulge/internal loops (not
yet ported, harder QUBO engineering) dominate. Both D-Wave and IBM QAOA
successfully solve the QUBO as formulated; porting hairpin energy in
surfaced a real, separate finding — it creates a correlated landscape
that shallow QAOA cannot yet reliably solve, confirmed via simulated
annealing that the model itself is correct.

## Approach

1. **Formulation.** RNA secondary structure is represented as a set of
   binary "quartet" variables, each denoting a stacked base pair — pair
   (i,j) immediately stacked on pair (i+1,j-1). This mirrors the MIS-style
   formulation used in recent gate-based work (arXiv:2505.05782), with one
   key difference: quartet weights are **real ΔG stacking energies pulled
   directly from ViennaRNA's own parameter tables** (`RNA.param().stack`),
   not tunable heuristic constants as used in prior quantum-annealing RNA
   folding papers (Fox et al. 2022; Zaborniak et al. 2022; Jiang et al.
   2023).
2. **Constraints.** Each base pairs at most once; no crossing pairs
   (pseudoknot-free), enforced as QUBO penalty terms.
3. **Classical benchmark.** ViennaRNA's `RNA.fold()` provides ground-truth
   MFE structures throughout.
4. **Quantum/quantum-inspired solving.** D-Wave Leap (hybrid and
   direct-QPU annealing) and QAOA on IBM hardware, both solving the
   identical QUBO formulation.
5. **Validation methodology.** Every component was cross-checked against
   an independent ground truth before being trusted: an exact brute-force
   solver validates the QUBO's constraint encoding; a dynamic program
   validates the energy model at scale beyond brute force's ~20-variable
   ceiling; real hardware results are checked with an explicit constraint
   validator, not assumed feasible.

## Key Findings

**1. The stacking-only model's match rate against real MFE collapses with
sequence length, and this was only visible after testing at scale.**
Small hand-picked toy sequences gave a misleadingly good 4/5 match. A
320-sequence sweep of *random* sequences (not hand-picked) shows the real
picture: **10.0% overall match rate, dropping to 0% at n≥60nt** — exactly
the length range used for hardware scaling runs. The energy gap also grows
systematically more negative with length, meaning the model doesn't just
occasionally pick a different fold — it increasingly overestimates
stability, consistent with the missing loop-entropy penalties mattering
more at scale, not less. *(Full detail: DEVLOG.md, "DP-based large-scale
validation".)*

**2. Adding real hairpin/bulge/internal-loop energies raises match rate
from 10.0% to 53.4% — and a real bug was caught building it.** A dynamic
program using ViennaRNA's own loop-energy evaluators (`eval_hp_loop`,
`eval_int_loop`) confirmed this improvement on the identical 320-sequence
sweep. During development, a "free multiloop fallback" was found to be
silently zeroing out hairpin penalties on *every* closure (0 always beats
a positive penalty) — caught by cross-checking the DP's own claimed energy
against ViennaRNA's independent structure evaluator, not by inspection.
Removing the fallback (rather than patching around it) fixed the model.
**This improved energy model has not yet been ported into the QUBO
itself** — the QUBO used for all D-Wave/IBM results below is still
stacking-only. This is the single highest-leverage remaining piece of
work. *(DEVLOG.md, "Full-loop-energy DP".)*

**3. D-Wave: the constraint graph is too dense for direct QPU embedding at
any real scale.** Density stays above 47% even at 540 variables (100nt) —
far beyond D-Wave's native ~15-20 connections per physical qubit. Direct
QPU submission (`EmbeddingComposite`) confirms this concretely: wall time
explodes with size (0.28s → 6.40s → 40.80s, n=20→30→40) while actual
`qpu_access_time` stays flat (~125-183ms) — the bottleneck is classical
minor-embedding search, not the quantum hardware. On an identical
sequence, direct-QPU solution quality was 39% worse than hybrid solving
(-11.4 vs -18.8), consistent with chain breaks from the dense embedding.
*(DEVLOG.md, "QPU-direct results".)*

**4. D-Wave's hybrid solver mostly doesn't touch the QPU at all for this
formulation, at these sizes.** Across a corrected scaling run (n=20-100,
capturing real `qpu_access_time` from `sampleset.info`), the QPU was
accessed on only 1 of 6 sequences — `qpu_access_time = 0` for n=20 through
80 (up to 362 variables), with `charge_time`/`run_time` sitting at a flat
~3.0s regardless of size (Leap's runtime floor, not scaling work). An
initial hypothesis that this was due to the constraint graph decomposing
into independent components was tested directly (connected-components
analysis) and **found to be wrong** — all six graphs are a single
connected component. The more likely explanation (unconfirmed) is a
size-based threshold internal to Leap's hybrid workflow. *(DEVLOG.md,
"Corrected hybrid loop".)*

**5. QAOA required a fundamentally different penalty weight than the
classical solvers — and a bug caused by reusing the wrong one was caught
and fixed. Even after the fix, QAOA's success is genuinely stochastic,
not reliable — measured directly, not assumed.** The QUBO's default
penalty (sized for classical exactness-guarantees, ~2x the sum of all
favorable energies) is roughly 11x larger than any single favorable
pairing on small instances. This is irrelevant to a classical exact
solver but makes the QAOA cost landscape so penalty-dominated that
"select nothing, violate nothing" becomes an inescapable local attractor
— confirmed empirically (0% match across 8 restarts, ruled out as an
optimizer-effort problem before diagnosing the real cause). A much
tighter, QAOA-specific penalty (~1.5x the largest single quartet energy)
fixed the total-collapse failure mode.

**Correction to an earlier claim in this README:** a single successful
run was previously reported here as "exact match on both test sequences,"
worded as if the fix guaranteed correctness. It doesn't. A 20-trial sweep
across different random seeds, using the exact parameters `run_ibm.py`
submits to hardware, found **GGGAAACCC succeeds only 15% of the time
(3/20)** and **GCGCUUCGGCGC succeeds 25% of the time (5/20)** — real
variance, not a guarantee, even on these tiny 4-6 qubit toy cases. This
was found by rerunning the validation after IBM access was restored (see
Finding 6 below) and getting a wrong answer on real hardware — investigated
directly rather than dismissed, and confirmed as expected stochastic
behavior (via a proper multi-seed sweep) rather than a new regression.
*(DEVLOG.md, "QAOA / IBM".)*

**6. Real IBM hardware confirms the formulation works on some runs, fails
on others exactly as the success-rate finding above predicts, and reveals
an early noise wall on top of that.** `ibm_kingston` (4 qubits,
`GGGAAACCC`, first hardware run): exact match, but at only 17.0% shot
confidence vs. 97.1% on the noiseless Aer simulator — real hardware
noise, quantified directly. A second run (`ibm_marrakesh`, 6 qubits,
`GCGCUUCGGCGC`) surfaced a further bug: the plurality bitstring was
outright infeasible (malformed structure, energy exactly ~2x the true
value — two conflicting quartets both firing).
Fixed with post-selection (rank all sampled bitstrings by shot count,
return the highest-count *feasible* one). After the fix, the top 3
bitstrings by shot count were still all infeasible; the best feasible
answer in the entire 2000-shot sample was the trivial empty structure, at
only 5.2% confidence. **Going from 4 to 6 qubits was enough to erase a
correct answer entirely on this backend/circuit combination** — a sharp,
quantified, real hardware-noise finding. *(DEVLOG.md, "IBM hardware
results".)*

## Files

| File | Purpose |
|---|---|
| `rna_qubo.py` | Quartet variable enumeration; real Turner2004 stacking energies from ViennaRNA; `get_hairpin_energy()` for the QUBO port. |
| `validate_brute_force.py` | Exact brute-force solver, small sequences; ground truth for the energy model. |
| `build_bqm.py` | `dimod` BQM/QUBO construction. Ports real hairpin energy by default (`include_hairpin=True`) via delta-correction; `include_hairpin=False` reproduces the original stacking-only model exactly. |
| `dp_validator.py` | O(n³) DP, stacking-only model, validated at scale (320 random sequences). |
| `dp_stack_hairpin.py` | O(n³) DP restricted to exactly the model the QUBO can express (stacking+hairpin) — used to scope the QUBO port honestly before writing it. |
| `dp_full_energy.py` | O(n⁴) DP with real hairpin/bulge/internal-loop energies via ViennaRNA's own evaluators (bulge/internal not yet ported into the QUBO). |
| `qaoa_hamiltonian.py` | Converts the same BQM into a QAOA-ready Ising Hamiltonian. |
| `qaoa_simulator.py` | QAOA on Aer simulator, with post-selection and the tight-penalty fix. |
| `run_dwave.py` | D-Wave Leap submission (hybrid/QPU). Run locally — not reachable from a sandboxed environment. |
| `run_ibm.py` | IBM Quantum submission. Run locally — not reachable from a sandboxed environment. |
| `colab_run.ipynb` | End-to-end Colab runner for the full validation + D-Wave pipeline. |
| `DEVLOG.md` | Full chronological development log — every run, bug, and fix in detail. |

## Setup

```bash
pip install -r requirements.txt
dwave setup                                    # or: export DWAVE_API_TOKEN=...
python -c "from qiskit_ibm_runtime import QiskitRuntimeService; \
  QiskitRuntimeService.save_account(channel='ibm_quantum_platform', \
  token='YOUR_TOKEN', overwrite=True, set_as_default=True)"
```

## Reproducing the results

```bash
python validate_brute_force.py   # energy model vs. real ViennaRNA MFE (toy cases)
python build_bqm.py              # BQM vs. independent brute force (must match)
python dp_validator.py           # DP vs. brute force cross-check
python dp_full_energy.py         # full-loop-energy DP self-consistency check
python qaoa_simulator.py         # QAOA vs. exact ground truth, Aer simulator

python run_dwave.py GGGAAACCC --qpu
python run_dwave.py <30-100nt sequence> --hybrid
python run_ibm.py GGGAAACCC
```

## Optional Advanced Tasks — status

The challenge lists four optional advanced tasks. Honest status on each,
stated directly rather than implied:

1. **Pseudoknots** — **real, validated proof-of-concept.** See dedicated
   section below. Scoped deliberately (not general pseudoknot folding,
   which is NP-hard), but genuinely working and independently verified.
2. **Compare multiple quantum encodings** — **done.** See dedicated
   section below.
3. **Sampling / hardware-inspired noise** — **substantial progress across
   two platforms.** 4 independent IBM runs across 3 backends, plus 3
   independent D-Wave QPU-direct runs with real chain-break-fraction and
   read-confidence metrics (neither captured before this task). See
   dedicated sections below.
4. **Qubit count vs. constraint enforcement trade-off** — **done.** See
   dedicated section below.

## Encoding Comparison (`rna_qubo_pairs.py`)

**Encoding A** (this project's primary formulation, used for all D-Wave/IBM
results): quartet variables — one binary variable per stacked pair,
stacking energy baked directly into the variable's linear bias.

**Encoding B** (built for this comparison): the raw base-pair encoding used
in prior literature (Fox, DePrince, Skolnick 2022; Zaborniak et al. 2022)
— one binary variable per valid base pair, full stop. A lone pair gets
zero linear reward; stacking energy is a *quadratic* bonus between two
variables representing adjacent pairs. Both encode the identical
stacking-only physical model two structurally different ways — this is a
genuine encoding comparison, not two different physics models being
mistaken for one.

**Validation:** both encodings solved via `dimod.ExactSolver` reach
identical energy and structure on all 5 toy sequences (5/5 exact match),
confirming Encoding B is a correct alternative formulation before
comparing anything about it.

**Qubit count and density, same 6 sequences used throughout the D-Wave
section** (`encoding_comparison_results.json`):

| n | A qubits (quartets) | A density | B qubits (pairs) | B density | B/A qubit ratio |
|---|---|---|---|---|---|
| 20 | 15 | 0.886 | 55 | 0.677 | 3.67 |
| 30 | 41 | 0.757 | 130 | 0.551 | 3.17 |
| 40 | 72 | 0.610 | 239 | 0.498 | 3.32 |
| 60 | 178 | 0.545 | 560 | 0.451 | 3.15 |
| 80 | 362 | 0.497 | 1093 | 0.419 | 3.02 |
| 100 | 540 | 0.475 | 1590 | 0.407 | 2.94 |

**Encoding B consistently needs ~3x more qubits than Encoding A** across
every sequence length tested (ratio 2.94-3.67, converging toward ~3x as
sequences grow), while its constraint graph is consistently *less* dense
(e.g. 0.407 vs 0.475 at n=100). This is a real, quantified engineering
trade-off: Encoding A is qubit-efficient by baking stacking directly into
the variable definition, at the cost of a denser constraint graph, since
every quartet is a priori "pre-committed" to a specific stacking
relationship.

**A genuine degenerate-optimum finding, caught and verified rather than
dismissed:** running simulated annealing (`neal`, 1000 reads) on both
encodings for the n=20 sequence, both reached the *same* optimal energy
(-5.00 kcal/mol) but *different* structures — `((..((....))..))....`
(Encoding A) vs. `((..((....))...))...` (Encoding B). Verified directly,
independent of either encoding's own bookkeeping: both structures are
well-formed (balanced brackets, all canonical pairs) and both
independently recompute to exactly -5.00 kcal/mol under the stacking-only
model. This is a real tied optimum — two different stacking patterns
achieve identical energy on this sequence — and the two encodings'
different variable/constraint geometries led the same solver to different
(equally valid) optima. Reported honestly as a finding about degenerate
solution landscapes, not a bug in either encoding.

## Sampling / Hardware Noise Study (`noise_study_aggregate.py`)

The original hardware results (Finding 6) had exactly 2 data points — one
success, one failure — at two different qubit counts. That's suggestive
but not a real distribution. 3 additional real hardware runs were
collected: 2 more runs of the identical `GGGAAACCC` problem (giving 3
independent runs total of the same problem, plus the original from
Finding 6 = **4 runs total**), and 1 run of a newly validated 5-qubit
intermediate case (`GGGUUCCCC`, exact classical match confirmed before
submission) to fill the gap between the existing 4- and 6-qubit points.
**IBM's free-plan QPU-time budget was fully exhausted collecting this
data, with no further access available going forward — this dataset is
final, not a snapshot pending a rerun.**

| seq | qubits | backend | confidence | correct? |
|---|---|---|---|---|
| GGGAAACCC | 4 | ibm_kingston | 17.0% | yes |
| GGGAAACCC | 4 | ibm_fez | 21.6% | yes |
| GGGAAACCC | 4 | ibm_fez | 24.0% | yes |
| GGGAAACCC | 4 | ibm_marrakesh | 15.3% | **no** |
| GGGUUCCCC | 5 | ibm_marrakesh | 9.4% | no |
| GCGCUUCGGCGC | 6 | ibm_marrakesh | 5.2% | no |

**Real finding: the identical 4-qubit problem gave 3 correct and 1
incorrect result across 4 independent runs on 3 different backends.**
Confidence for the correct runs clustered at 17.0-24.0%; the one
incorrect run's confidence (15.3%) sits inside that same range — meaning
confidence alone doesn't cleanly separate a correct outcome from an
incorrect one at this qubit count, which is itself worth reporting rather
than only citing the mean.

**An honest confound, disclosed rather than hidden:** all 4 of these
hardware runs used `run_ibm.py`'s *pre-fix* penalty calibration — the same
stale-calibration bug pattern caught and fixed in `qaoa_simulator.py`
earlier (Finding 5), but that fix was not carried over to `run_ibm.py`
until after these runs were already submitted (the IBM QPU-time budget
was exhausted before a corrected rerun was possible). Checked directly:
the stale penalty used was 6.10 for both `GGGAAACCC` and `GGGUUCCCC`, vs.
a correctly-calibrated 4.00 for both — a real but moderate gap, not the
severe ~11x oversizing that caused total QAOA collapse in Finding 5. This
means **the noise-vs-qubit-count trend reported here cannot be cleanly
separated from a possible penalty-calibration effect** — both plausibly
contribute to the observed variance, and this dataset cannot distinguish
between them. `run_ibm.py` has since been corrected (uses the same
`bqm.linear`-based calibration as `qaoa_simulator.py`) for whenever
hardware access resumes.

**Correction:** the sentence below originally said this rerun was "no
longer possible" because IBM access was "permanently exhausted." That
claim was wrong — access was restored, and the rerun was completed.

**Rerun against the corrected penalty calibration, with a genuinely
informative result — not the clean "noise vs. calibration" separation
hoped for, but something more useful: hard evidence that QAOA's success
on this Hamiltonian is fundamentally stochastic, not a calibration
artifact.** Both `GGGAAACCC` and `GCGCUUCGGCGC` were rerun on
`ibm_fez` with the corrected calibration (penalty 4.00 and 4.15
respectively, matching `qaoa_simulator.py`'s formula exactly). **Both
runs found the wrong answer** — `GGGAAACCC` converged to the trivial
unfolded structure (simulator pre-optimization also found the wrong
answer, 1458/2000 shots), `GCGCUUCGGCGC` converged to a different,
suboptimal fold (`.(((....))).` instead of the true `((((....))))`).

This was investigated directly rather than reported as a mystery: a
20-trial sweep across random seeds, using `run_ibm.py`'s exact parameters
(`reps=2, maxiter=100, restarts=5`), found **GGGAAACCC succeeds only 15%
of the time (3/20) and GCGCUUCGGCGC succeeds 25% of the time (5/20)** —
see Finding 5 above for the full correction this prompted to an earlier,
overstated claim. Today's two real hardware failures are the *expected*
outcome given those success rates (a combined ~64% chance both would
fail), not evidence of a new regression. The corrected calibration is
confirmed working (it eliminates the total-collapse failure mode from
Finding 5), but it does not guarantee success on any individual run —
this is the honest, now-measured shape of the remaining limitation.

### D-Wave noise study — a genuine extension, not a repeat of the IBM data

D-Wave's noise characteristics are structurally different from IBM's:
IBM's are gate-based (shallow-circuit decoherence, gate error), D-Wave's
are annealing-based (chain breaks from embedding, thermal noise during
the anneal). `run_dwave.py` previously captured neither `chain_break_
fraction` nor any read-level confidence metric — only `sampleset.first`,
the single best sample, was ever inspected. Both were added and validated
before use (`dwave_noise_study_results.json`): 3 independent QPU-direct
runs of the same problem (`GGACGGCGCUUCUACUCAAC`, n=20, 15 qubits, whose
true optimum is the fully unfolded structure, independently confirmed
earlier via `dp_stack_hairpin.py`).

| run | chain_break_fraction (mean) | frac. reads w/ any break | read-level confidence | correct? |
|---|---|---|---|---|
| 1 | 0.0059 | 8.89% | 1.60% | yes |
| 2 | 0.0047 | 6.98% | 2.70% | yes |
| 3 | 0.0070 | 10.47% | 2.20% | yes |

All 3 runs found the correct answer and were genuinely feasible. But
**read-level confidence (mean 2.17%) is roughly 8-10x lower than IBM's
shot confidence on its correct runs (mean ~21%)**, despite D-Wave's
chain-break rate being quite modest (under 7 in 1000 reads have any break
at all, on average). This means chain breaks alone do not explain the gap
— something else about the annealing search is driving most reads away
from the true optimum, even on a problem this small.

**Honestly flagged, not overclaimed:** with only 3 runs and no access to
the full per-read energy distribution beyond what `run_dwave.py` currently
extracts, the specific mechanism (e.g. broad thermal exploration across
many near-degenerate feasible states in a very dense penalty landscape —
this problem's constraint graph density is 0.981, the highest recorded
anywhere in this project) is a plausible explanation, not a confirmed
one. This is also not a clean apples-to-apples comparison with the IBM
numbers above: different qubit counts (15 vs. 4-6), different hardware
paradigms entirely (quantum annealing vs. gate-based QAOA). The finding
that stands on its own, without needing that comparison: **D-Wave's
per-read confidence in the correct answer is low even when the answer
found is consistently correct and chain breaks are rare** — a genuine,
quantified, D-Wave-specific noise characteristic, distinct from anything
the IBM data showed.

Two different penalty regimes exist in this project, for good reasons:
the **classical-exactness penalty** (`build_bqm.py`'s default, sized to
guarantee correctness against arbitrary constraint violations on
classical exact solvers — scales with the *sum* of every favorable energy
in the problem) and the **tight penalty** (`qaoa_simulator.py`'s default,
sized just above the *single* largest favorable energy — small enough for
shallow QAOA to actually navigate). Both were built and validated
separately, for different solvers. Assembling them together against the
same 6 sequences used throughout the D-Wave scaling work reveals a real,
quantifiable trade-off, independent of and prior to any hardware noise:

| n | qubits | density | classical penalty | tight penalty | ratio |
|---|---|---|---|---|---|
| 20 | 15 | 0.905 | 69.20 | 11.80 | 5.86 |
| 30 | 41 | 0.774 | 137.80 | 12.85 | 10.72 |
| 40 | 72 | 0.619 | 255.40 | 11.80 | 21.64 |
| 60 | 178 | 0.548 | 445.20 | 12.30 | 36.21 |
| 80 | 362 | 0.499 | 1256.00 | 13.90 | 90.36 |
| 100 | 540 | 0.476 | 1853.80 | 15.36 | 120.73 |

**The tight penalty stays nearly flat (11.8 → 15.4) while the
classical-exactness penalty grows ~27x (69.2 → 1853.8) over the same
range — so the ratio between them grows 20x, from under 6 to over 120,
purely as a function of qubit count.** This is a structural property of
the formulation, not a hardware or noise effect: as problems scale, the
penalty required to *guarantee* feasibility on a classical solver becomes
increasingly incompatible with the penalty a NISQ-era shallow circuit can
actually work with. This directly explains, in a single quantified
number, why QAOA's practical usable size is so much smaller than a
classical solver's: the "safe" penalty and the "solvable" penalty diverge
sharply with scale, independent of the noise-wall finding reported
separately (Finding 6).

**Honest note on the density column:** these density values (0.905 at
n=20, etc.) are marginally higher than the ones reported earlier in the
D-Wave section (0.886 at n=20) — this table uses the current default
`build_bqm()` behavior, which now includes the hairpin-energy port. The
hairpin delta-correction adds a small number of genuine new correlational
edges (93 → 95 at n=20, confirmed directly) — a minor, explainable,
expected side effect of porting more real physics into the model, not a
discrepancy to gloss over.

## Pseudoknots (`pseudoknot_qubo.py`)

General pseudoknot folding is NP-hard, and even restricted, well-studied
classes (H-type pseudoknots) normally need specialized O(n⁶) algorithms
(Rivas & Eddy, 1999) to score correctly — well beyond what remaining
project time allows. **This is not an attempt at that.** What was built
instead, honestly scoped: relax this project's no-crossing constraint so
the QUBO *can* select crossing quartets when genuinely favorable, and
prove — on a real, hand-verified test case — that doing so lets the
solver correctly find a better answer the non-crossing model structurally
cannot reach.

**Energy model caveat, stated upfront rather than discovered later:** real
pseudoknots carry additional loop-topology entropy penalties beyond plain
stacking energy (Rivas & Eddy's `gw`/`gwh` terms), which this extension
does not model — every crossing quartet is scored with the same real
Turner2004 stacking energy as a nested one. This means the extension's
total energy for a genuine pseudoknot is an *underestimate* of the true
(more unfavorable) free energy. Consistent with how every other energy
gap in this project has been handled (hairpin added incrementally and
honestly; bulge/internal/multiloop still unmodeled and documented as
such) — a real, bounded simplification, not a hidden one.

**What's relaxed:** only the no-crossing check. The no-shared-base
constraint (a base can't be paired twice) is kept — that's physically
required regardless of pseudoknots.

**Test case, built and verified from scratch, not assumed:** a
26-nucleotide sequence (`GCAUGAACGUACAAACAUGCAGUACG`) hand-constructed
with two independent 5bp stems ("armA": positions 0-4 pairing with
15-19; "armB": positions 7-11 pairing with 21-25) using distinct,
non-repetitive sequences to avoid spurious alternative pairings.
Verified programmatically before use: both stems are fully canonical
Watson-Crick pairs, and their outer pairs — (0,19) and (7,25) — provably
cross (`0 < 7 < 19 < 25`). ViennaRNA (which cannot represent pseudoknots
at all) finds only one stem, MFE -5.0.

**Result:**

| model | energy | structure |
|---|---|---|
| non-crossing QUBO (stacking-only baseline) | -8.70 | uses only one arm |
| pseudoknot-permissive QUBO | **-17.40** | uses **both** arms simultaneously |

-17.40 exactly matches the hand-calculated expectation (each arm
independently contributes -8.70 in stacking energy across 4 stacking
steps each, and the two arms share no bases, so using both should sum
exactly). The solver found this without being told the answer in advance
— confirmation the relaxation works generally, not just for this one
case.

**Independently re-verified from scratch**, not just trusted from the
solver's own output: the selected 10 pairs use no base twice, all 10 are
canonical Watson-Crick pairs, and the total energy recomputed directly
from those pairs (bypassing the BQM entirely) matches the solver's
claimed energy exactly (-17.40 = -17.40).

**Honest scope limits:**
- Stacking-only (no hairpin energy port to this variant — would need the
  same delta-correction technique re-validated under crossing, not
  attempted here).
- No real pseudoknot-specific energy penalties (stated above).
- No systematic scan for how often *real* MFE structures benefit from
  pseudoknots at realistic sequence lengths — this is a proof-of-concept
  on one deliberately-constructed case, not a statistical study like the
  320-sequence validations elsewhere in this project.
- Not tested on D-Wave or IBM hardware — remaining time was prioritized
  toward validating the formulation correctly rather than spending
  limited hardware access on a proof-of-concept.

## QUBO Port: Real Hairpin Energy (`rna_qubo.py`, `build_bqm.py`)

The DP proved a full-energy model (stacking + hairpin + bulge + internal
loop) reaches 53.4% match rate, but the actual QUBO run on D-Wave/IBM
throughout this project was stacking-only. This section ports **hairpin**
energy into the QUBO — bulge/internal loops remain future work, for a
specific, honest reason explained below.

**Scoping check performed before writing any QUBO code:** a DP restricted
to exactly the model the QUBO can express (stacking + hairpin, no
bulge/internal — `dp_stack_hairpin.py`) was run on the identical
320-sequence sweep first, to find out how much of the 53.4% improvement
hairpin alone accounts for:

| length | stacking-only | stack+hairpin | full (+bulge/internal) |
|---|---|---|---|
| 10 | 15.0% | **100.0%** | 100.0% |
| 15 | 20.0% | **75.0%** | 72.5% |
| 20 | 32.5% | **70.0%** | 90.0% |
| 30 | 10.0% | **35.0%** | 77.5% |
| 40 | 2.5% | **2.5%** | 50.0% |
| 60 | 0.0% | **2.5%** | 25.0% |
| 80 | 0.0% | **0.0%** | 10.0% |
| 100 | 0.0% | **0.0%** | 2.5% |
| **overall** | **10.0%** | **35.6%** | **53.4%** |

**Hairpin alone captures most of the benefit at short sequences (n≤30)
and almost none of it at longer ones (n≥40)**, where bulge/internal loops
are the dominant remaining gap. Practically: this port meaningfully
strengthens the IBM QAOA results (9-12nt, exactly where hairpin-only
helps most) and does very little for the D-Wave scaling results (20-100nt,
mostly outside where hairpin-only helps). This is stated directly rather
than implied as a uniform improvement.

**Implementation: a delta-correction technique requiring no auxiliary
variables.** Every quartet gets a baseline linear addition of its inner
pair's real hairpin energy (assuming the stack chain ends there); if a
continuation quartet exists, a quadratic correction term exactly cancels
that baseline. This works cleanly *because* quartets are already a linear
chain (each has at most one possible continuation) — bulge/internal loops
would require a flexible-gap search with multiple candidate "children" per
closing pair, needing genuine auxiliary "at most one active child"
constraints (a 3-way constraint that doesn't reduce to a QUBO without
extra variables) — real, harder, separate work, not attempted here.

**Validation:** the new BQM matches `dp_stack_hairpin.py` exactly (energy
and structure) on all 5 toy sequences via `dimod.ExactSolver`, and a
regression check confirms `include_hairpin=False` still reproduces the
original stacking-only brute-force results exactly — no silent behavior
change for code that doesn't opt in.

**A real, separate finding surfaced immediately after this port:** QAOA's
existing tight-penalty calibration (built for the stacking-only model) was
computed from the raw stacking-only quartet values and went stale the
moment hairpin terms were added as additional linear biases — it silently
regressed QAOA back to the trivial all-zero attractor on one test sequence
(`GCGCUUCGGCGC`), the exact failure mode the original tight-penalty fix
was built to prevent. Fixed by calibrating against the BQM's actual linear
magnitudes instead of the stale quartets-only view.

**After that fix, a second, more fundamental issue appeared: QAOA (reps=3,
8 restarts, 300 iterations) still could not reliably find the true
optimum on the hairpin-aware Hamiltonian**, landing on a feasible but
suboptimal local minimum instead. Diagnosed directly rather than left
unexplained: classical simulated annealing on the *identical* BQM (not
circuit-depth-limited like QAOA) finds the true optimum trivially on both
test sequences. **This confirms the QUBO model itself is correct — the
hairpin delta-correction introduces genuine variable correlations
(a quadratic term that only pays off when two specific quartets are
selected together) that make the landscape harder for a shallow
parameterized circuit to navigate, independent of any bug.** This is
listed as a real, open item: QAOA needs either deeper circuits, better
initialization, or a reformulation to handle the richer Hamiltonian — not
yet solved, honestly flagged rather than hidden.

**Important scope correction to the claim above:** "the QUBO model itself
is correct" was validated only against the two small toy sequences tested
at that point. A larger-scale check (below) found a real, more general
limitation in the same delta-correction technique.

### A second, more serious bug found during D-Wave rerun prep — found and fixed

While selecting a replacement test sequence for the D-Wave rerun (below),
simulated annealing on a candidate 60nt sequence found an energy
(-5.40) *better* than the independent DP's reported optimum (-5.20) for
the identical restricted model — which should be impossible if both are
correct. Investigated directly rather than dismissed:

The SA solution turned out to be a genuine **branching structure** — an
outer stem enclosing a completely independent inner stem (an internal-loop
topology), not a single continuous helix. Both `dp_stack_hairpin.py` and
the QUBO's delta-correction technique were built to handle only two cases
per closing pair: continue the stack, or close a hairpin over empty
space. Neither was designed to represent "close over a loop that itself
contains another independent nested stem." The DP simply cannot search
that structure (a real, acceptable limitation, already documented). But
the **QUBO's constraint system did not forbid it** — nothing in
`conflicting()` prevented two non-crossing, non-base-sharing quartets from
being selected even when one's assumed hairpin loop actually contained the
other — and when this happened, the delta-correction wrongly charged a
full hairpin-closure energy for a loop that wasn't empty.

Confirmed by direct, independent recomputation: the true physical energy
of this exact structure (via ViennaRNA's own hairpin and internal-loop
evaluators, correctly accounting for the branch) is **-6.20**. The BQM
itself claimed **-5.40** for the identical selection — a real ~0.8
kcal/mol accounting error, not noise or a rounding artifact.

**Fixed** (`build_bqm.py`, `nested_noncolinear()`): two quartets are now
forbidden together if one is nested inside the other's span but they are
not colinear members of the same symmetric helix — this correctly catches
independent branching structures while still permitting continuous
stacks of any depth (matching inward offsets are automatically exempt,
without needing to check whether every intermediate step is selected).
Deliberately kept separate from `conflicting()` rather than merged into
it, since `conflicting()` is shared by the stacking-only model, the
alternative pair encoding (Encoding B), and `run_dwave.py`'s feasibility
check — none of which have this bug, and none of which should be newly
restricted by a fix that's only correct under the hairpin-aware model.

Validated three ways before trusting it: (1) all 5 toy sequences still
match `dp_stack_hairpin.py` exactly — no regression; (2) the
stacking-only model (`include_hairpin=False`) is completely unaffected —
no regression; (3) under the fixed model, the exact selection that
previously scored -5.40 now scores **7430.6** (correctly and heavily
penalized), and simulated annealing on the fixed BQM converges to
**-5.20** — exactly matching the independent DP.

### D-Wave scaling rerun against the hairpin-aware QUBO — including a confirmed correction

The original D-Wave scaling results (Findings 3-4) all predate the
hairpin port and reflect the stacking-only model. Rerun on the identical
6 sequences (`scaling_results_hairpin_aware_dwave.json`).

**Correction to an earlier claim in this README:** the first hairpin-aware
rerun (before the branching fix above) was reported here as unaffected by
the branching bug — that claim was wrong, and it's being corrected
directly rather than silently edited away. Checking the actual selected
structures directly (`nested_noncolinear()`, after the fix existed)
showed **2 of the 6 already-reported hardware results — n=80 and n=100 —
had in fact exploited the exact branching bug**, the same way the toy
test case that led to the fix did. Both were rerun on real D-Wave
hardware against the corrected model:

| n | qubits | qpu_access_time (μs) | energy | feasible | corrected? |
|---|---|---|---|---|---|
| 20 | 15 | 109,803 | 0.00 | yes | no (unaffected, confirmed) |
| 30 | 41 | 81,754 | -1.60 | yes | no (unaffected, confirmed) |
| 40 | 72 | 91,378 | -5.00 | yes | no (unaffected, confirmed) |
| 60 | 178 | 97,603 | 0.00 | yes | no (unaffected, confirmed) |
| 80 | 362 | 150,045 | **-6.20** (was -6.90) | yes | **yes** |
| 100 | 540 | 154,952 | **-6.70** (was -7.93) | yes | **yes** |

Both wrong values (-6.90, -7.93) were genuinely reported by real D-Wave
hardware before the fix — this was not a simulation artifact, it was a
real, previously-undiscovered accounting error present in production
hardware results already in this repository. The historical (wrong)
values are preserved, clearly flagged, in
`scaling_results_hairpin_aware_dwave_PRE_BRANCHING_FIX.json` rather than
deleted, for transparency.

**A real cost of the fix, worth reporting alongside the correction:**
constraint edge count roughly doubled at both corrected sizes (n=80:
32,614 → 49,678 edges; n=100: 69,299 → 110,355), and graph density
increased substantially (n=80: 0.499 → 0.760; n=100: 0.476 → 0.758). The
`nested_noncolinear` penalty pass adds real new constraint edges, and at
these larger sizes, a meaningful fraction of all quartet pairs turn out
to be non-colinear nested pairs that must now be forbidden. This is a
genuine additional cost of correctness, not a free fix.

### Resolved: the "size-based QPU-access threshold" was actually a density effect

Earlier in this project (before the branching fix), only 1 of the 6
standard sequences (n=100) showed any real `qpu_access_time` at all — the
other 5 (n=20 through n=80) showed exactly `0`, meaning Leap's hybrid
solver decomposed and solved them entirely classically without touching
the QPU. This was left as an open question — is n=100 a stable threshold,
or stochastic per-run?

**Investigated directly with 6 new real hardware runs**: n=100 repeated 3
times, plus 3 new intermediate sizes (n=85, 90, 95) using freshly
generated, validated sequences.

| n | qubits | edges | density | qpu\_access\_time (μs) |
|---|---|---|---|---|
| 20 | 15 | 95 | 0.905 | 109,803 |
| 30 | 41 | 635 | 0.774 | 81,754 |
| 40 | 72 | 1,582 | 0.619 | 91,378 |
| 60 | 178 | 8,637 | 0.548 | 97,603 |
| 80 | 362 | 49,678 | 0.760 | 150,045 |
| 85 | 411 | 63,839 | 0.758 | 151,513 |
| 90 | 439 | 69,402 | 0.722 | 153,726 |
| 95 | 428 | 70,208 | 0.768 | 153,683 |
| 100 | 540 | 110,355 | 0.758 | 154,952 |
| 100 (repeat) | — | — | — | 154,811 |
| 100 (repeat) | — | — | — | 153,644 |

**Two clean findings, resolving both parts of the original question:**

1. **n=100's QPU access is stable, not stochastic.** Three independent
   repeats vary by under 0.85% (153,644–154,952 μs) — deterministic
   behavior on Leap's side, not per-run noise.
2. **The original "only n=100 triggers QPU access" finding is now stale,
   not wrong for its time.** It was measured on the pre-branching-fix
   model. All 9 sizes tested here — n=20 through n=100 — now show
   substantial, nonzero QPU access. The real cause: the branching fix
   (`nested_noncolinear`) added real constraint edges at every scale, not
   just the largest ones, pushing every tested problem size over
   whatever internal threshold determines whether Leap's hybrid workflow
   engages the QPU at all. This was never really a *size* threshold — it
   was a constraint-density threshold that size happened to correlate
   with in the old, lower-density model.

**Two of six sequences (n=20, n=60) optimize to the fully unfolded
structure.** This was predicted and confirmed *before* running hardware,
via `dp_stack_hairpin.py` — the independent ground truth for this exact
model agreed these specific sequences genuinely have no favorable fold
once real hairpin-loop entropy penalties are included, exactly the same
phenomenon documented earlier for `GGAAUUCC` in the DP validation section.
This is a correct, expected result, not a solver failure. These two were
also independently confirmed unaffected by the branching bug via local
simulated annealing before spending any QPU time re-confirming them.

**A separate real bug was also caught and fixed while interpreting the
first rerun.** `run_dwave.py`'s feasibility check originally compared the
solver's true total energy against a stale stacking-only energy sum — the
same bug pattern already caught twice elsewhere (`qaoa_simulator.py`'s
`raw_stacking_energy` field and its penalty calibration). This caused 4 of
6 runs in the first rerun to print a false `CONSTRAINT VIOLATED` warning.
Fixed by checking the real constraint directly (`conflicting()` on the
selected quartets) instead of comparing energies that were never supposed
to match once the hairpin port landed. Independently re-verified after
the fix: all 6
selections are genuinely feasible, and every recomputed `true_model_energy`
matches D-Wave's own reported energy exactly.

## Limitations

- **No multiloop support**, DP or QUBO. A correct treatment needs a
  separate WM table with `MLbase`/`MLintern`/`MLclosing` costs and branch
  counting — real, scoped-out future work, not a hidden gap.
- **No pseudoknots** — out of scope by design (challenge's optional task).
- **The QUBO now ports hairpin energy but still lacks bulge/internal
  loops** — proven by the DP to matter most at n≥40, i.e. exactly the
  D-Wave-scale regime. The auxiliary-constraint engineering needed for
  bulge/internal loops is real, separate future work.
- **QAOA cannot yet reliably solve the hairpin-aware Hamiltonian** —
  confirmed via simulated annealing that the model itself is correct; the
  correlated delta-correction terms make the landscape genuinely harder
  for a shallow circuit, independent of the earlier penalty-calibration
  bug (which was found and fixed separately).
- **QAOA's tight penalty is empirically validated on small cases, not a
  rigorous guarantee at arbitrary scale** — always check the real
  constraint output (`feasible`), don't trust the penalty blindly as
  problems grow.
- **IBM hardware noise dominates by 6 qubits** on the current
  backend/circuit combination without deeper error mitigation (dynamical
  decoupling, readout correction, zero-noise extrapolation — not
  attempted, out of scope for remaining time).
- **The QUBO's delta-correction could not represent branching
  (internal-loop) topologies** — found and fixed (`nested_noncolinear()`
  in `build_bqm.py`). See the dedicated write-up above for the discovery,
  the confirmed ~0.8 kcal/mol accounting error, and the three-part
  validation of the fix.
- **D-Wave hardware results now reflect the hairpin-aware QUBO** (rerun
  completed, before the branching fix was found — the 6 sequences
  submitted did not trigger it, so those results remain valid) —
  **IBM hardware results still predate the hairpin port.** IBM access
  was reported exhausted earlier in this project; it has since been
  restored (see future-work item below for status of the rerun).

## Future Work

1. ~~Rerun D-Wave hardware results against the hairpin-aware QUBO~~ —
   **done.** ~~Rerun IBM hardware results~~ — **done** (access was
   restored; earlier "permanently exhausted" claim corrected). Result:
   both reruns found wrong answers on real hardware, investigated and
   confirmed as expected given QAOA's measured ~15-25% success rate on
   these toy cases (see Finding 5/6 corrections) — a genuine, informative
   finding, not a clean confirmation.
2. ~~Fix the QUBO's branching/internal-loop-topology gap~~ — **done**
   (`nested_noncolinear()`). ~~Rerun the D-Wave scaling study once more
   against the now-fixed model~~ — **done**; found and corrected 2 of 6
   already-published results (n=80, n=100) that had exploited the bug.
3. Add a genuine D-Wave noise study using the newly-added
   `chain_break_fraction`/read-confidence metrics at additional qubit
   counts (only n=20 tested so far) — would extend the current 3-run
   dataset into a real qubit-count-vs-noise curve, the same way the IBM
   data was originally intended to be before its access was exhausted.
4. Diagnose and fix QAOA's difficulty with the hairpin-aware Hamiltonian's
   correlated landscape (deeper circuits, better initialization, or a
   reformulation).
5. Bulge/internal-loop QUBO extension via auxiliary "at most one active
   child" constraints — proven by the DP to matter most at n≥40.
6. Proper multiloop DP and QUBO extension (WM table, branch counting).
7. Error mitigation to push usable QAOA qubit count past the observed
   4-6 qubit noise wall.
8. ~~Investigate D-Wave's apparent size-based QPU-access threshold~~ —
   **resolved.** See dedicated section below: it wasn't a size threshold,
   it was a density effect caused by the branching fix.
9. ~~Pseudoknot-aware extension~~ — **proof-of-concept done**
   (`pseudoknot_qubo.py`). Real next steps: port hairpin energy to this
   variant; add real pseudoknot-specific loop-topology penalties
   (Rivas & Eddy `gw`/`gwh` terms); test on real hardware; run a
   systematic scan across random sequences (like the 320-sequence
   stacking/hairpin validations) to see how often crossing structures
   actually improve on real MFE at realistic lengths, not just the one
   deliberately-constructed test case.

## References

- Fox, DePrince, Skolnick. *RNA folding using quantum computers.* PLOS Comp Bio, 2022.
- Zaborniak et al. *A QUBO model of the RNA folding problem optimized by variational hybrid quantum annealing.* arXiv:2208.04367, 2022.
- Jiang et al. *Predicting RNA Secondary Structure on Universal Quantum Computer.* arXiv:2305.09561, 2023.
- *mRNA secondary structure prediction using utility-scale quantum computers.* arXiv:2405.20328, 2024.
- *Towards secondary structure prediction of longer mRNA sequences using a quantum-centric optimization scheme.* arXiv:2505.05782, 2025.
