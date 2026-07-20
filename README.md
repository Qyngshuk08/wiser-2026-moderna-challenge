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
10% of the time on unbiased random sequences (0% at n≥60nt). Adding real
hairpin/bulge/internal-loop energies (validated via a classical dynamic
program, not yet ported into the QUBO itself) raises this to 53.4%. Both
D-Wave and IBM QAOA successfully solve the QUBO as formulated; the
remaining gap to real MFE is a property of the energy model, not the
solver.

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
| `rna_qubo.py` | Quartet variable enumeration; real Turner2004 stacking energies from ViennaRNA. |
| `validate_brute_force.py` | Exact brute-force solver, small sequences; ground truth for the energy model. |
| `build_bqm.py` | `dimod` BQM/QUBO construction with constraint penalties; cross-validated against brute force. |
| `dp_validator.py` | O(n³) DP, stacking-only model, validated at scale (320 random sequences). |
| `dp_full_energy.py` | O(n⁴) DP with real hairpin/bulge/internal-loop energies via ViennaRNA's own evaluators. |
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

## Limitations

- **No multiloop support**, DP or QUBO. A correct treatment needs a
  separate WM table with `MLbase`/`MLintern`/`MLclosing` costs and branch
  counting — real, scoped-out future work, not a hidden gap.
- **No pseudoknots** — out of scope by design (challenge's optional task).
- **The QUBO itself is stacking-only** even though the DP proves a
  full-energy model is achievable and worth porting — this is the single
  biggest gap between what's validated and what's running on hardware.
- **QAOA's tight penalty is empirically validated on small cases, not a
  rigorous guarantee at arbitrary scale** — always check the real
  constraint output (`feasible`), don't trust the penalty blindly as
  problems grow.
- **IBM hardware noise dominates by 6 qubits** on the current
  backend/circuit combination without deeper error mitigation (dynamical
  decoupling, readout correction, zero-noise extrapolation — not
  attempted, out of scope for remaining time).

## Future Work

1. Port hairpin/bulge/internal-loop energies into the QUBO itself
   (highest priority — proven worthwhile by the DP, not yet done).
2. Proper multiloop DP and QUBO extension.
3. Error mitigation techniques to push the usable QAOA qubit count past
   the observed 4-6 qubit noise wall.
4. Investigate the D-Wave hybrid solver's apparent size-based QPU-access
   threshold directly (repeat n=100, test intermediate sizes).
5. Pseudoknot-aware extension (challenge's optional advanced task).

## References

- Fox, DePrince, Skolnick. *RNA folding using quantum computers.* PLOS Comp Bio, 2022.
- Zaborniak et al. *A QUBO model of the RNA folding problem optimized by variational hybrid quantum annealing.* arXiv:2208.04367, 2022.
- Jiang et al. *Predicting RNA Secondary Structure on Universal Quantum Computer.* arXiv:2305.09561, 2023.
- *mRNA secondary structure prediction using utility-scale quantum computers.* arXiv:2405.20328, 2024.
- *Towards secondary structure prediction of longer mRNA sequences using a quantum-centric optimization scheme.* arXiv:2505.05782, 2025.
