"""
Optional advanced task: "Evaluate the approach under sampling or
hardware-inspired noise."

Aggregates every IBM hardware run recorded so far (ibm_hardware_results_
run*.json) into a single noise-vs-qubit-count table. Where multiple runs
exist for the same (seq, backend) pair, reports mean/std of shot
confidence instead of a single point -- a real study needs a distribution,
not two anecdotes, which is exactly what the original 2-point result
(Findings 6 in the README) was honest about not being.

Run this AFTER collecting additional real hardware runs (see
NOISE_STUDY_PLAN in this file's __main__ block for the exact commands) --
it only aggregates existing JSON result files, it does not touch hardware
itself (not reachable from this sandbox, same limitation as run_ibm.py).
"""

import glob
import json
import statistics


def load_all_results():
    records = []
    for path in sorted(glob.glob("ibm_hardware_results_run*.json")):
        with open(path) as f:
            data = json.load(f)
            for r in data:
                r["_source_file"] = path
                records.append(r)
    return records


def aggregate(records):
    # Deduplicate exact-duplicate records first. Real hardware noise means
    # two genuinely independent runs essentially never produce identical
    # shot counts -- an exact numeric match (same seq, backend, shots,
    # total_shots) is almost certainly the same execution recorded twice
    # for context, not a second real run. Confirmed directly: run2.json's
    # GGGAAACCC entry duplicates run1.json's, added as a comparison
    # reference when GCGCUUCGGCGC was recorded -- not an independent run.
    seen = set()
    deduped = []
    for r in records:
        shots = r.get("shots") or r.get("chosen_shots") or r.get("best_bitstring_shots")
        key = (r["seq"], r.get("backend"), shots, r.get("total_shots"))
        if key in seen:
            continue
        seen.add(key)
        deduped.append(r)

    grouped = {}
    for r in deduped:
        key = (r["seq"], r.get("backend", "unknown"))
        grouped.setdefault(key, []).append(r)

    summary = []
    for (seq, backend), runs in grouped.items():
        confidences = []
        for r in runs:
            shots = r.get("shots") or r.get("chosen_shots") or r.get("best_bitstring_shots")
            total = r.get("total_shots", 2000)
            if shots is not None:
                confidences.append(shots / total)
        n_qubits = runs[0].get("num_qubits")
        matched = [r.get("match_exact", r.get("feasible")) for r in runs]
        summary.append({
            "seq": seq, "n": len(seq), "backend": backend,
            "num_qubits": n_qubits, "num_runs": len(runs),
            "mean_confidence": round(statistics.mean(confidences), 4) if confidences else None,
            "stdev_confidence": round(statistics.stdev(confidences), 4) if len(confidences) > 1 else None,
            "confidences": [round(c, 4) for c in confidences],
        })
    return sorted(summary, key=lambda x: x["num_qubits"] or 0)


def aggregate_by_sequence(records):
    """Group by sequence ONLY, ignoring backend -- for this noise study,
    different backends on the same logical problem is itself part of the
    finding (e.g. GGGAAACCC ran on both ibm_fez and ibm_marrakesh with
    materially different outcomes), and splitting them into separate rows
    (as aggregate() does, grouped by (seq, backend)) would hide that."""
    seen = set()
    deduped = []
    for r in records:
        shots = r.get("shots") or r.get("chosen_shots") or r.get("best_bitstring_shots")
        key = (r["seq"], r.get("backend"), shots, r.get("total_shots"))
        if key in seen:
            continue
        seen.add(key)
        deduped.append(r)

    grouped = {}
    for r in deduped:
        grouped.setdefault(r["seq"], []).append(r)

    summary = []
    for seq, runs in grouped.items():
        confidences = []
        per_run = []
        for r in runs:
            shots = r.get("shots") or r.get("chosen_shots") or r.get("best_bitstring_shots")
            total = r.get("total_shots", 2000)
            match = r.get("match_exact", r.get("feasible"))
            backend = r.get("backend", "unknown")
            if shots is not None:
                conf = shots / total
                confidences.append(conf)
                per_run.append({"backend": backend, "confidence": round(conf, 4), "match_exact": match})
        summary.append({
            "seq": seq, "n": len(seq), "num_qubits": runs[0].get("num_qubits"),
            "num_runs": len(runs),
            "mean_confidence": round(statistics.mean(confidences), 4) if confidences else None,
            "stdev_confidence": round(statistics.stdev(confidences), 4) if len(confidences) > 1 else None,
            "min_confidence": round(min(confidences), 4) if confidences else None,
            "max_confidence": round(max(confidences), 4) if confidences else None,
            "per_run": per_run,
        })
    return sorted(summary, key=lambda x: x["num_qubits"] or 0)


if __name__ == "__main__":
    records = load_all_results()
    if not records:
        print("No ibm_hardware_results_run*.json files found in this directory.")
    else:
        seq_summary = aggregate_by_sequence(records)
        print(f"{'seq':<15}{'qubits':>8}{'n_runs':>8}{'mean_conf':>12}{'min':>8}{'max':>8}")
        for s in seq_summary:
            print(f"{s['seq']:<15}{s['num_qubits']:>8}{s['num_runs']:>8}"
                  f"{s['mean_confidence']:>12.4f}{s['min_confidence']:>8.4f}{s['max_confidence']:>8.4f}")
            for pr in s["per_run"]:
                print(f"    -> {pr['backend']:<16} confidence={pr['confidence']:.4f}  match_exact={pr['match_exact']}")
        with open("noise_study_summary.json", "w") as f:
            json.dump(seq_summary, f, indent=2)

    print()
    print("=" * 70)
    print("NOISE STUDY PLAN -- commands to run in Colab with IBM hardware access")
    print("=" * 70)
    print("""
Currently have (from earlier sessions):
  - GGGAAACCC     (4 qubits): 1 run,  confidence 17.0%
  - GCGCUUCGGCGC  (6 qubits): 1 run,  confidence  5.2% (post-selection fallback)

To get a real study instead of 2 anecdotes, run:

  # 3 more runs on the existing 4-qubit case, for a real distribution
  !python run_ibm.py GGGAAACCC
  !python run_ibm.py GGGAAACCC
  !python run_ibm.py GGGAAACCC

  # NEW: validated 5-qubit intermediate point (fills the 4->6 gap)
  !python run_ibm.py GGGUUCCCC
  !python run_ibm.py GGGUUCCCC
  !python run_ibm.py GGGUUCCCC

After each run, append its JSON output to a new
ibm_hardware_results_run3.json (same schema as run1/run2), then run:

  !python noise_study_aggregate.py

to get the assembled mean/std-by-qubit-count table.
""")
