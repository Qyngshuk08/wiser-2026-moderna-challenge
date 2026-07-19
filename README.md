# WISER 2026 — Moderna Challenge: Quantum RNA Secondary Structure Prediction

Team: **Qudit Creons**

Optimization of mRNA secondary structure prediction (MFE folding) using a
quantum/quantum-inspired QUBO formulation, benchmarked against ViennaRNA.

## Approach

Prior quantum-annealing RNA folding work (Fox et al. 2022, PLOS Comp Bio;
Zaborniak et al. 2022, arXiv:2208.04367; Jiang et al. 2023, arXiv:2305.09561;
2025 quartet formulation, arXiv:2505.05782) formulates RNA folding as
maximizing base-pair or stem count, using **tunable heuristic reward
constants**, not real thermodynamic parameters.

This project instead pulls **actual Turner2004 nearest-neighbor stacking
energies directly from ViennaRNA's own parameter tables**
(`RNA.param().stack`) and uses those as QUBO weights. The result is a
formulation whose objective is a real (partial) free energy, not an
arbitrary proxy for one.

**Documented limitation:** this first-pass model captures stacking energy
only. Hairpin/bulge/internal-loop/multiloop entropic penalties are
length-dependent and not naturally quadratic, so they are not yet included.
See `validate_brute_force.py` output for a concrete, quantified example of
where this causes a structure mismatch against true MFE, and by how much
the omitted loop penalty matters.

## Files

| File | Purpose |
|---|---|
| `rna_qubo.py` | Enumerates valid base pairs and stacked-pair ("quartet") variables; pulls real Turner2004 stacking energies from ViennaRNA. |
| `validate_brute_force.py` | Exact brute-force solver over feasible quartet subsets (one-pair-per-base, no pseudoknots); compares against ViennaRNA's true MFE structure. Ground-truth check on the energy model itself, independent of any QUBO/solver machinery. |
| `build_bqm.py` | Builds the actual `dimod` BQM/QUBO with penalty terms enforcing the same constraints; cross-validated against `validate_brute_force.py` (5/5 exact match on toy sequences). |
| `run_dwave.py` | D-Wave Leap submission script (hybrid or QPU). Must be run locally with a valid `DWAVE_API_TOKEN` — not runnable in a sandboxed/offline environment. |

## Setup

```bash
pip install ViennaRNA dimod dwave-ocean-sdk
dwave setup      # or: export DWAVE_API_TOKEN=your_token
dwave ping       # confirm connectivity before submitting real jobs
```

## Validate the model (no hardware required)

```bash
python validate_brute_force.py   # energy model vs. real ViennaRNA MFE
python build_bqm.py              # BQM vs. independent brute force (must match)
```

Current result: 4/5 toy sequences match ViennaRNA's MFE dot-bracket exactly.
The one mismatch (`GGAAUUCC`) is a weak single-stack stem that our model
folds but real MFE leaves unfolded — the hairpin-loop entropy penalty (not
yet modeled) outweighs the stacking bonus for marginal stems. This is
reported as a known, quantified boundary of the current formulation, not
hidden.

## Run on D-Wave

```bash
python run_dwave.py GGGAAACCC --qpu                 # sanity check, ~4 variables
python run_dwave.py <30-60nt sequence> --hybrid      # scaling study
```

Each run reports variable count, constraint-edge count, wall time, energy,
resulting structure, and whether the constraint penalty was satisfied
(feasible solution). These numbers are the basis for the
scaling/quantum-resource analysis deliverable.

## Status / Next steps

- [x] Real-energy stacking-only QUBO, validated against ViennaRNA on toy sequences
- [x] BQM construction validated against independent brute force
- [ ] D-Wave Leap hybrid scaling runs across sequence lengths
- [ ] IBM QPU QAOA comparison at small scale
- [ ] Dynamic-programming-scale validation (beyond brute-force's ~20-variable limit)
- [ ] Optional: loop-length-binned energy terms (hairpin/bulge/internal loop) as an extension

## References

- Fox, DePrince, Skolnick. *RNA folding using quantum computers.* PLOS Comp Bio, 2022.
- Zaborniak et al. *A QUBO model of the RNA folding problem optimized by variational hybrid quantum annealing.* arXiv:2208.04367, 2022.
- Jiang et al. *Predicting RNA Secondary Structure on Universal Quantum Computer.* arXiv:2305.09561, 2023.
- *mRNA secondary structure prediction using utility-scale quantum computers.* arXiv:2405.20328, 2024.
- *Towards secondary structure prediction of longer mRNA sequences using a quantum-centric optimization scheme.* arXiv:2505.05782, 2025.
