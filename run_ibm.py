"""
Run this LOCALLY, not in the sandbox that built it -- IBM Quantum's runtime
API isn't reachable from here.

Setup (once):
    pip install qiskit-ibm-runtime
    # save your IBM Quantum API token -- IBM retired the old "ibm_quantum"
    # channel; only ibm_cloud and ibm_quantum_platform are valid now:
    python -c "from qiskit_ibm_runtime import QiskitRuntimeService; \
        QiskitRuntimeService.save_account(channel='ibm_quantum_platform', \
        token='YOUR_TOKEN', overwrite=True, set_as_default=True)"
    # if that alone errors, you likely also need your instance CRN from the
    # 'Instances' tab on your IBM Quantum Platform dashboard:
    #   QiskitRuntimeService.save_account(channel='ibm_quantum_platform',
    #       token='YOUR_TOKEN', instance='YOUR_CRN', overwrite=True,
    #       set_as_default=True)

Usage:
    python run_ibm.py GGGAAACCC
    python run_ibm.py GCGCUUCGGCGC --backend ibm_fez

Small sequences only -- realistically 10-15 qubits before circuit depth and
noise make QAOA output unusable on real hardware. Uses the tight-penalty
fix from qaoa_simulator.py (penalty ~1.5x max single quartet energy, NOT
the classical-exactness penalty from build_bqm.py) or QAOA will very likely
collapse to the trivial all-zero solution on real hardware too, same as it
did on the simulator before the penalty was fixed.

NOTE on timing: job.metrics()['usage']['quantum_seconds'] is captured on
the FINAL 2000-shot sampling job only, not on each of the ~50 COBYLA
optimizer iterations beforehand (each of those also submits a real,
separate job -- adding per-iteration metrics capture would add overhead
for a number that's largely redundant with the final job's). This is the
real QPU execution time for one representative job at this problem size,
not the total across every job the run submits. Neither prior hardware
run (run1, run2) captured this at all -- don't assume or estimate a
number for those; they simply don't have one.
"""

import argparse
import sys
sys.path.insert(0, ".")

from qiskit_ibm_runtime import QiskitRuntimeService, SamplerV2 as RuntimeSamplerV2
from qiskit.circuit.library import QAOAAnsatz
from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager
from scipy.optimize import minimize
import numpy as np

from qaoa_hamiltonian import build_qaoa_hamiltonian
from validate_brute_force import pairs_used, to_dot_bracket
from build_bqm import conflicting
from itertools import combinations


def run(seq, backend_name=None, reps=2, maxiter=50, penalty=None):
    if penalty is None:
        preview = build_qaoa_hamiltonian(seq, penalty=1.0)["quartets"]
        max_single = max(abs(e) for e in preview.values()) if preview else 1.0
        penalty = 1.5 * max_single + 1.0
        print(f"using tight penalty={penalty:.2f} (NOT the classical-exactness "
              f"penalty -- see qaoa_simulator.py notes on why)")

    result = build_qaoa_hamiltonian(seq, penalty=penalty)
    op = result["ising_op"]
    offset = result["ising_offset"]
    var_list = result["var_list"]
    quartets = result["quartets"]

    print(f"seq length: {len(seq)}   qubits needed: {result['num_qubits']}")

    service = QiskitRuntimeService(channel="ibm_quantum_platform")
    backend = service.backend(backend_name) if backend_name else service.least_busy(operational=True, simulator=False)
    print(f"backend: {backend.name}")

    ansatz = QAOAAnsatz(cost_operator=op, reps=reps)
    ansatz.measure_all()

    pm = generate_preset_pass_manager(optimization_level=1, backend=backend)
    isa_ansatz = pm.run(ansatz)
    isa_op = op.apply_layout(isa_ansatz.layout)

    sampler = RuntimeSamplerV2(mode=backend)

    x0 = np.random.default_rng(42).uniform(0, 2 * np.pi, ansatz.num_parameters)

    # NOTE: real hardware queue time makes a full COBYLA optimization loop
    # expensive (each iteration = one queued job). This does a SHORT
    # optimization budget by default -- raise maxiter only if you have
    # queue time to spare, and expect this to take a while either way.
    def cost_fn(params):
        bound = isa_ansatz.assign_parameters(params)
        job = sampler.run([bound], shots=1000)
        result = job.result()[0]
        counts = result.data.meas.get_counts()
        # estimate expectation from counts against the ORIGINAL (unmapped) op
        total = sum(counts.values())
        # cheap proxy: use the best bitstring's raw stacking energy as the
        # objective for the classical optimizer, rather than a full
        # expectation value recompute per iteration (expensive on hardware)
        best_bits = max(counts, key=counts.get)[::-1]
        assignment = {v: int(best_bits[i]) for i, v in enumerate(var_list)}
        selected = [v for v, val in assignment.items() if val == 1]
        return sum(quartets[q] for q in selected)

    opt = minimize(cost_fn, x0, method="COBYLA", options={"maxiter": maxiter})

    bound = isa_ansatz.assign_parameters(opt.x)
    job = sampler.run([bound], shots=2000)
    counts = job.result()[0].data.meas.get_counts()

    # Real QPU execution time -- neither prior hardware run captured this,
    # and it should not be estimated or guessed. job.metrics() exposes the
    # actual usage.quantum_seconds IBM bills/reports for this job.
    try:
        metrics = job.metrics()
        quantum_seconds = metrics.get("usage", {}).get("quantum_seconds")
    except Exception as e:
        metrics = None
        quantum_seconds = None
        print(f"WARNING: could not retrieve job.metrics(): {e}")

    print(f"job id: {job.job_id()}")
    print(f"real QPU execution time (quantum_seconds): {quantum_seconds}")

    # Post-select on feasibility instead of blindly taking the plurality
    # bitstring. Real hardware noise can make the top-shot answer
    # infeasible even with a tight penalty (confirmed empirically:
    # GCGCUUCGGCGC on ibm_marrakesh returned an infeasible top bitstring --
    # malformed dot-bracket, energy exactly ~2x the true value, consistent
    # with two conflicting quartets both firing). Rank all sampled
    # bitstrings by shot count and take the highest-count one that is
    # ACTUALLY feasible, rather than trusting the raw plurality blindly.
    ranked = sorted(counts.items(), key=lambda kv: -kv[1])
    chosen = None
    for bitstring, shots in ranked:
        bits = bitstring[::-1]
        assignment = {v: int(bits[i]) for i, v in enumerate(var_list)}
        selected = [v for v, val in assignment.items() if val == 1]
        is_feasible = not any(conflicting(q1, q2) for q1, q2 in combinations(selected, 2))
        if is_feasible:
            chosen = (bitstring, shots, selected, is_feasible)
            break
    if chosen is None:
        # no feasible bitstring in the entire sample -- report the raw
        # plurality anyway, but flag it loudly rather than hide it
        bitstring, shots = ranked[0]
        bits = bitstring[::-1]
        assignment = {v: int(bits[i]) for i, v in enumerate(var_list)}
        selected = [v for v, val in assignment.items() if val == 1]
        is_feasible = False
        chosen = (bitstring, shots, selected, is_feasible)
        print("WARNING: no feasible bitstring found anywhere in the 2000-shot "
              "sample. Reporting the raw plurality bitstring, which is "
              "infeasible. This is a real result, not a bug -- report it as-is.")

    best_bitstring, best_shots, selected, is_feasible = chosen
    rank_of_chosen = [b for b, _ in ranked].index(best_bitstring) + 1

    raw_energy = sum(quartets[q] for q in selected)
    db = to_dot_bracket(len(seq), pairs_used(selected))

    print(f"structure: {db}")
    print(f"raw stacking energy: {raw_energy:.2f}")
    print(f"feasible: {is_feasible}")
    print(f"chosen bitstring rank by shot count: {rank_of_chosen} "
          f"(1 = plurality; >1 means the plurality bitstring was infeasible "
          f"and this is the best feasible fallback)")
    print(f"chosen bitstring shots: {best_shots}/{sum(counts.values())}")

    return {
        "seq": seq, "n": len(seq), "backend": backend.name,
        "num_qubits": result["num_qubits"], "structure": db,
        "raw_energy": raw_energy, "feasible": is_feasible,
        "shots": best_shots, "total_shots": sum(counts.values()),
        "rank_by_shot_count": rank_of_chosen,
        "quantum_seconds": quantum_seconds, "job_id": job.job_id(),
    }


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("sequence")
    ap.add_argument("--backend", default=None)
    ap.add_argument("--reps", type=int, default=2)
    ap.add_argument("--maxiter", type=int, default=50)
    args = ap.parse_args()
    run(args.sequence.upper(), args.backend, args.reps, args.maxiter)
