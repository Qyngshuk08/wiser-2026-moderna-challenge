"""
Run this LOCALLY, not in the sandbox that built it -- D-Wave's cloud API
isn't reachable from here.

Setup (once):
    pip install dwave-ocean-sdk
    dwave setup            # or: export DWAVE_API_TOKEN=your_leap_token

Usage:
    python run_dwave.py GGGGAAAACCCCUUUUGGGGAAAACCCC --hybrid
    python run_dwave.py GGGGAAAACCCCUUUUGGGGAAAACCCC --qpu

Use --hybrid (LeapHybridSampler) for anything above ~30-40nt -- this is
where your scaling study should live, since hybrid solvers absorb much
larger QUBOs than the QPU alone.

Use --qpu (EmbeddingComposite + DWaveSampler) only for small sequences
where you want a genuine QPU-only result to report, not hybrid-solved.
Expect embedding to fail or need heavy chain-strength tuning above ~15-20
variables on a single QPU -- that failure point IS scaling data, report it
even if it's ugly.
"""

import argparse
import time
import sys
sys.path.insert(0, ".")  # rna_qubo.py / build_bqm.py must be alongside this file

from build_bqm import build_bqm
from validate_brute_force import pairs_used, to_dot_bracket


def run(seq, mode, min_loop=3, num_reads=1000):
    bqm, quartets = build_bqm(seq, min_loop)
    print(f"seq length: {len(seq)}   quartet variables: {len(bqm.variables)}   "
          f"quadratic terms (constraint edges): {len(bqm.quadratic)}")

    if mode == "hybrid":
        from dwave.system import LeapHybridSampler
        sampler = LeapHybridSampler()
        t0 = time.time()
        sampleset = sampler.sample(bqm)
        wall_time = time.time() - t0
    else:
        from dwave.system import DWaveSampler, EmbeddingComposite
        sampler = EmbeddingComposite(DWaveSampler())
        t0 = time.time()
        sampleset = sampler.sample(bqm, num_reads=num_reads)
        wall_time = time.time() - t0

    best = sampleset.first
    selected = [q for q, v in best.sample.items() if v == 1]
    db = to_dot_bracket(len(seq), pairs_used(selected))

    print(f"mode: {mode}   wall time: {wall_time:.2f}s")
    print(f"best energy (incl. penalty, should be near pure-stacking energy "
          f"if feasible): {best.energy:.2f}")
    print(f"structure: {db}")

    raw_energy = sum(quartets[q] for q in selected)
    print(f"raw stacking energy of selected quartets: {raw_energy:.2f}  "
          f"({'feasible' if abs(raw_energy - best.energy) < 1e-6 else 'CONSTRAINT VIOLATED -- raise penalty weight'})")

    return {
        "seq": seq, "n": len(seq), "mode": mode,
        "num_variables": len(bqm.variables), "num_edges": len(bqm.quadratic),
        "wall_time_s": wall_time, "energy": best.energy,
        "structure": db, "feasible": abs(raw_energy - best.energy) < 1e-6,
    }


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("sequence")
    ap.add_argument("--hybrid", action="store_true")
    ap.add_argument("--qpu", action="store_true")
    ap.add_argument("--min-loop", type=int, default=3)
    ap.add_argument("--num-reads", type=int, default=1000)
    args = ap.parse_args()
    mode = "hybrid" if args.hybrid else "qpu"
    run(args.sequence.upper(), mode, args.min_loop, args.num_reads)
