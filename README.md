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
why — including two real bugs caught mid-project by cross-checking results
against independent ground truth rather than trusting plausible-looking
output.

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
and fixed.** The QUBO's default penalty (sized for classical
exactness-guarantees, ~2x the sum of all favorable energies) is roughly
11x larger than any single favorable pairing on small instances. This is
irrelevant to a classical exact solver but makes the QAOA cost landscape
so penalty-dominated that "select nothing, violate nothing" becomes an
inescapable local attractor — confirmed empirically (0% match across 8
restarts, ruled out as an optimizer-effort problem before diagnosing the
real cause). A much tighter, QAOA-specific penalty (~1.5x the largest
single quartet energy) fixed this: exact match on both test sequences,
confirmed with a real constraint check, not an energy-gap proxy.
*(DEVLOG.md, "QAOA / IBM".)*

**6. Real IBM hardware confirms the formulation works, and reveals an
early noise wall.** `ibm_kingston` (4 qubits, `GGGAAACCC`): exact match,
but at only 17.0% shot confidence vs. 97.1% on the noiseless Aer
simulator — real hardware noise, quantified directly. A second run
(`ibm_marrakesh`, 6 qubits, `GCGCUUCGGCGC`) surfaced a further bug: the
plurality bitstring was outright infeasible (malformed structure, energy
exactly ~2x the true value — two conflicting quartets both firing).
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

1. **Pseudoknots** — not attempted. Excluded by design throughout (the
   no-crossing constraint is structural, not incidental).
2. **Compare multiple quantum encodings** — not attempted. All hardware
   results use one encoding (quartets/stacked-pairs) on two backends
   (D-Wave, IBM); that is a backend comparison, not an encoding
   comparison.
3. **Sampling / hardware-inspired noise** — partial, not a designed study.
   Real IBM hardware results (see above) quantify actual noise (97.1% →
   17.0% → 5.2% shot confidence across two qubit counts), which is
   stronger evidence than a synthetic noise model, but was not a
   deliberately swept experiment.
4. **Qubit count vs. constraint enforcement trade-off** — partial, data
   exists but wasn't assembled as its own analysis: D-Wave's constraint
   graph density vs. embedding success (Finding 3) and QAOA's classical
   vs. tight penalty finding (Finding 5) are both directly relevant.

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
- **D-Wave/IBM hardware results in this repo predate the hairpin port** —
  all scaling and hardware runs reflect the original stacking-only QUBO.
  Rerunning them against the hairpin-aware model is real future work.

## Future Work

1. Rerun D-Wave/IBM hardware results against the hairpin-aware QUBO —
   biggest expected impact at the short sequence lengths used in the IBM
   QAOA runs.
2. Diagnose and fix QAOA's difficulty with the hairpin-aware Hamiltonian's
   correlated landscape (deeper circuits, better initialization, or a
   reformulation).
3. Bulge/internal-loop QUBO extension via auxiliary "at most one active
   child" constraints — proven by the DP to matter most at n≥40.
4. Proper multiloop DP and QUBO extension (WM table, branch counting).
5. Error mitigation to push usable QAOA qubit count past the observed
   4-6 qubit noise wall.
6. Investigate D-Wave's apparent size-based QPU-access threshold directly
   (repeat n=100, test intermediate sizes).
7. Pseudoknot-aware extension (challenge's optional advanced task).
8. Build a second quantum encoding and compare against the quartet
   formulation (challenge's optional advanced task).
9. Assemble the existing qubit-count/constraint-enforcement data (D-Wave
   density findings + QAOA penalty findings) into a dedicated trade-off
   analysis (challenge's optional advanced task).

## References

- Fox, DePrince, Skolnick. *RNA folding using quantum computers.* PLOS Comp Bio, 2022.
- Zaborniak et al. *A QUBO model of the RNA folding problem optimized by variational hybrid quantum annealing.* arXiv:2208.04367, 2022.
- Jiang et al. *Predicting RNA Secondary Structure on Universal Quantum Computer.* arXiv:2305.09561, 2023.
- *mRNA secondary structure prediction using utility-scale quantum computers.* arXiv:2405.20328, 2024.
- *Towards secondary structure prediction of longer mRNA sequences using a quantum-centric optimization scheme.* arXiv:2505.05782, 2025.
