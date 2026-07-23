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
from itertools import combinations
sys.path.insert(0, ".")  # rna_qubo.py / build_bqm.py must be alongside this file

from build_bqm import build_bqm, conflicting
from validate_brute_force import pairs_used, to_dot_bracket


def run(seq, mode, min_loop=3, num_reads=1000):
    bqm, quartets = build_bqm(seq, min_loop)
    v = len(bqm.variables)
    e = len(bqm.quadratic)
    max_possible = v * (v - 1) // 2
    density = e / max_possible if max_possible else 0.0
    print(f"seq length: {len(seq)}   quartet variables: {v}   "
          f"quadratic terms (constraint edges): {e}   "
          f"graph density: {density:.3f}  "
          f"(D-Wave QPU native connectivity is ~15-20 per qubit -- dense "
          f"logical graphs like this embed poorly or not at all without "
          f"the hybrid solver's classical decomposition)")

    if mode == "hybrid":
        from dwave.system import LeapHybridSampler
        sampler = LeapHybridSampler()
        t0 = time.time()
        sampleset = sampler.sample(bqm)
        wall_time = time.time() - t0
        chain_break_stats = None
    else:
        from dwave.system import DWaveSampler, EmbeddingComposite
        sampler = EmbeddingComposite(DWaveSampler())
        t0 = time.time()
        sampleset = sampler.sample(bqm, num_reads=num_reads)
        wall_time = time.time() - t0

        # Real noise metrics for QPU-direct mode -- neither was captured
        # before. chain_break_fraction is specific to embedded QPU solves
        # (irrelevant to hybrid, which doesn't need chain embedding) and
        # is a genuine, D-Wave-specific noise signal, structurally
        # different from IBM's gate-based shot-confidence noise.
        try:
            cbf = sampleset.record.chain_break_fraction
            chain_break_stats = {
                "mean": float(cbf.mean()), "max": float(cbf.max()),
                "min": float(cbf.min()), "frac_with_any_break": float((cbf > 0).mean()),
            }
        except Exception as e:
            chain_break_stats = None
            print(f"WARNING: could not extract chain_break_fraction: {e}")

    best = sampleset.first
    selected = [q for q, v in best.sample.items() if v == 1]
    db = to_dot_bracket(len(seq), pairs_used(selected))

    # Read-level confidence: what fraction of all reads landed on the
    # BEST energy found, vs scattered across worse answers? This is the
    # D-Wave analog of IBM's shot-confidence metric -- not previously
    # computed here at all (only .first, the single best sample, was
    # ever inspected).
    read_confidence = None
    if mode == "qpu":
        try:
            total_reads = int(sampleset.record.num_occurrences.sum())
            best_energy_reads = int(sum(
                rec.num_occurrences for rec in sampleset.record
                if abs(rec.energy - best.energy) < 1e-6
            ))
            read_confidence = best_energy_reads / total_reads if total_reads else None
        except Exception as e:
            print(f"WARNING: could not compute read-level confidence: {e}")

    # sampleset.info carries the actual quantum-resource breakdown -- wall
    # time alone is dominated by classical presolve on hybrid solves and is
    # NOT a scaling signal by itself. Field names differ slightly between
    # LeapHybridSampler and QPU-direct solves, so pull whatever's present.
    info = dict(sampleset.info) if hasattr(sampleset, "info") else {}
    qpu_access_time_us = info.get("qpu_access_time")  # hybrid: microseconds
    charge_time_us = info.get("charge_time")
    run_time_us = info.get("run_time")
    timing = info.get("timing", {})  # QPU-direct solves nest timing here

    print(f"mode: {mode}   wall time: {wall_time:.2f}s")
    print(f"sampleset.info: {info}")
    print(f"best energy (incl. penalty, should be near pure-stacking energy "
          f"if feasible): {best.energy:.2f}")
    print(f"structure: {db}")
    if chain_break_stats:
        print(f"chain_break_fraction: mean={chain_break_stats['mean']:.4f}  "
              f"max={chain_break_stats['max']:.4f}  "
              f"frac_reads_with_any_break={chain_break_stats['frac_with_any_break']:.4f}")
    if read_confidence is not None:
        print(f"read-level confidence (fraction of {num_reads} reads at best energy): "
              f"{read_confidence:.4f}")

    # Real feasibility check, not a stale energy-gap proxy. The old check
    # compared against `sum(quartets[q] for q in selected)` -- the
    # STACKING-ONLY energy -- but best.energy now includes real hairpin
    # baseline/correction terms since build_bqm()'s include_hairpin=True
    # became the default. Those two were never going to match once the
    # hairpin port landed; this is the same stale-comparison bug pattern
    # already caught and fixed twice elsewhere (qaoa_simulator.py's
    # raw_stacking_energy, and its penalty calibration). Fixed the same
    # way: check the REAL constraint directly (no shared bases, no
    # crossing) rather than comparing energies that were never supposed
    # to be equal post-hairpin-port.
    is_feasible = not any(conflicting(q1, q2) for q1, q2 in combinations(selected, 2))
    true_model_energy = bqm.energy({v: (1 if v in selected else 0) for v in bqm.variables})
    raw_energy = sum(quartets[q] for q in selected)
    print(f"raw stacking-only energy of selected quartets (excludes hairpin terms): {raw_energy:.2f}")
    print(f"true total model energy (stacking + hairpin, via bqm.energy()): {true_model_energy:.2f}")
    print(f"feasible (real constraint check): {is_feasible}")

    def jsonable(x):
        # dimod/numpy scalars aren't JSON-serializable directly -- cast at
        # the source instead of patching it after every run.
        if isinstance(x, dict):
            return {str(k): jsonable(v) for k, v in x.items()}
        if hasattr(x, "item"):
            return x.item()
        return x

    return {
        "seq": seq, "n": len(seq), "mode": mode,
        "num_variables": len(bqm.variables), "num_edges": len(bqm.quadratic),
        "wall_time_s": jsonable(wall_time),
        "qpu_access_time_us": jsonable(qpu_access_time_us),
        "charge_time_us": jsonable(charge_time_us),
        "run_time_us": jsonable(run_time_us),
        "qpu_timing_raw": jsonable(timing),
        "energy": jsonable(best.energy),
        "structure": db, "feasible": bool(is_feasible),
        "raw_stacking_energy": raw_energy, "true_model_energy": true_model_energy,
        "chain_break_stats": chain_break_stats,
        "read_confidence": read_confidence,
        "num_reads": num_reads if mode == "qpu" else None,
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
