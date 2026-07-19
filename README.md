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
- [x] D-Wave Leap hybrid scaling runs, n=20 to n=100 (`scaling_results_run1.json`) — all 6 structures independently verified: balanced dot-bracket, all canonical pairs, min loop length respected
- [ ] IBM QPU QAOA comparison at small scale
- [ ] Dynamic-programming-scale validation (beyond brute-force's ~20-variable limit)
- [ ] Optional: loop-length-binned energy terms (hairpin/bulge/internal loop) as an extension

### QPU-direct results (`scaling_results_run2_qpu_direct.json`), n=20/30/40

Two findings, both stronger and more specific than "it doesn't scale":

3. **The bottleneck is classical embedding time, not the QPU.** Wall time: 0.28s → 6.40s → 40.80s (n=20→30→40). `qpu_access_time` over the same range: 126ms → 183ms → 125ms — flat, not growing. `EmbeddingComposite`'s minor-embedding search (minorminer) is almost certainly what's exploding, not the quantum hardware itself, and it's exploding because the constraint graph density (0.886 → 0.757 → 0.610 over this range) is far above the ~15-20 native connections per physical qubit that D-Wave's topology supports.
4. **Solution quality degrades on identical inputs.** The n=40 sequence here is the same sequence run under `--hybrid` earlier: hybrid found energy -11.4 [-18.8], a 39% worse result. Chain breaks from the dense embedding are the likely cause, though `chain_break_fraction` isn't currently captured by `run_dwave.py` — logged as a script improvement, not yet done.

**Outstanding:** the corrected hybrid loop (capturing `qpu_access_time` per size, not just wall time) has not yet been rerun after the Colab kernel caching issue was fixed — the six-sequence hybrid dataset in `scaling_results_run1.json` above still predates the `sampleset.info` capture fix and should eventually be redone for a clean side-by-side comparison, though the QPU-vs-hybrid quality gap on n=40 above is already sufficient evidence on its own.

### Corrected hybrid loop, n=20-100 (`scaling_results_run3_hybrid.json`) — the most important finding so far

5. **`qpu_access_time` is 0 for n=20, 30, 40, 60, and 80.** The QPU was not touched at all on 5 of 6 sequences, up to 362 variables. `charge_time`/`run_time` sit at ~2.99-3.0 million μs on every single run regardless of size — this is LeapHybridSampler's runtime floor, not work that scales with problem size. Only at n=100 does `qpu_access_time` become nonzero (103,707μs ≈ 104ms).

   The honest reading: **this formulation, at these sizes, is being solved entirely by the hybrid solver's classical presolve.** Looking at the resulting structures (mostly small, disjoint stem-loops rather than one large interacting fold) explains why — a problem that decomposes into small independent components is exactly what classical presolve handles without needing the QPU, regardless of how dense or how many total variables the QUBO has on paper. The n=100 case, where QPU access finally triggers, is the first (and so far only) genuinely interesting data point for actual quantum resource usage, and is worth investigating for what distinguishes it structurally from the rest — whether that's chance sequence composition or something about scale forcing less separable structure.

   This matters for how the whole project should be framed: raw variable/edge counts alone do not indicate quantum resource usage for this formulation. Whether the QPU gets invoked at all is a property of how separable the resulting structure is, not just problem size.

## References

- Fox, DePrince, Skolnick. *RNA folding using quantum computers.* PLOS Comp Bio, 2022.
- Zaborniak et al. *A QUBO model of the RNA folding problem optimized by variational hybrid quantum annealing.* arXiv:2208.04367, 2022.
- Jiang et al. *Predicting RNA Secondary Structure on Universal Quantum Computer.* arXiv:2305.09561, 2023.
- *mRNA secondary structure prediction using utility-scale quantum computers.* arXiv:2405.20328, 2024.
- *Towards secondary structure prediction of longer mRNA sequences using a quantum-centric optimization scheme.* arXiv:2505.05782, 2025.
