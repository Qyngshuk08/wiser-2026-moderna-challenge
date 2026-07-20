"""
Manual QAOA using QAOAAnsatz + Qiskit's current V2 primitives directly.
qiskit_algorithms 0.4.0 (latest available) still expects the V1 primitives
interface, which is incompatible with the installed Qiskit 2.5.0 -- rather
than pin an old Qiskit version and hope nothing else breaks, QAOA is
implemented directly here: it's a short enough algorithm that "roll it
ourselves against the still-working part of the API" is more reliable than
depending on an unmaintained compatibility shim.

Validated against dimod's ExactSolver ground truth before being trusted for
anything real-hardware-bound.
"""

import numpy as np
from scipy.optimize import minimize
from qiskit.circuit.library import QAOAAnsatz
from qiskit_aer.primitives import EstimatorV2, SamplerV2

from qaoa_hamiltonian import build_qaoa_hamiltonian
from build_bqm import solve_exact, conflicting
from validate_brute_force import pairs_used, to_dot_bracket
from itertools import combinations


def run_qaoa_simulator(seq, min_loop=3, reps=2, maxiter=100, seed=42, restarts=1, penalty=None):
    # QAOA needs a much tighter penalty than the classical exact solvers do.
    # build_bqm()'s default penalty (2x sum of ALL favorable energies) is
    # sized to guarantee exactness against ARBITRARY constraint violations
    # on ANY problem size -- but on a small QAOA instance it makes the
    # penalty ~10x larger than any single favorable energy, which traps
    # shallow QAOA at the trivial "violate nothing, gain nothing" all-zero
    # solution (confirmed empirically: default penalty -> 0% match across
    # 8 restarts on two test sequences; a tight penalty -> exact match on
    # both). Default here to just above the largest single quartet energy,
    # which is provably sufficient to forbid any single pairwise violation
    # while keeping the landscape far less rugged.
    if penalty is None:
        quartets_preview = build_qaoa_hamiltonian(seq, min_loop, penalty=1.0)["quartets"]
        max_single = max(abs(e) for e in quartets_preview.values()) if quartets_preview else 1.0
        penalty = 1.5 * max_single + 1.0

    result = build_qaoa_hamiltonian(seq, min_loop, penalty=penalty)
    op = result["ising_op"]
    offset = result["ising_offset"]
    var_list = result["var_list"]
    quartets = result["quartets"]

    cost_ansatz = QAOAAnsatz(cost_operator=op, reps=reps).decompose(reps=3)
    sample_ansatz = QAOAAnsatz(cost_operator=op, reps=reps).decompose(reps=3)
    sample_ansatz.measure_all()

    estimator = EstimatorV2()
    sampler = SamplerV2()

    def cost_fn(params):
        job = estimator.run([(cost_ansatz, op, params)])
        return job.result()[0].data.evs

    rng = np.random.default_rng(seed)
    best_opt = None
    for r in range(restarts):
        x0 = rng.uniform(0, 2 * np.pi, cost_ansatz.num_parameters)
        opt = minimize(cost_fn, x0, method="COBYLA", options={"maxiter": maxiter})
        if best_opt is None or opt.fun < best_opt.fun:
            best_opt = opt
    opt = best_opt

    bound_circuit = sample_ansatz.assign_parameters(opt.x)
    job = sampler.run([bound_circuit], shots=2000)
    counts = job.result()[0].data.meas.get_counts()

    # Post-select on feasibility, same fix as run_ibm.py after real
    # hardware demonstrated the plurality bitstring can be infeasible even
    # with a tight penalty.
    ranked = sorted(counts.items(), key=lambda kv: -kv[1])
    chosen = None
    for bitstring, shots in ranked:
        bits = bitstring[::-1]
        assignment = {v: int(bits[i]) for i, v in enumerate(var_list)}
        selected_try = [v for v, val in assignment.items() if val == 1]
        if not any(conflicting(q1, q2) for q1, q2 in combinations(selected_try, 2)):
            chosen = (bitstring, shots, selected_try)
            break
    if chosen is None:
        bitstring, shots = ranked[0]
        bits = bitstring[::-1]
        assignment = {v: int(bits[i]) for i, v in enumerate(var_list)}
        chosen = (bitstring, shots, [v for v, val in assignment.items() if val == 1])

    best_bitstring, best_shots, selected = chosen

    raw_energy = sum(quartets[q] for q in selected)
    db = to_dot_bracket(len(seq), pairs_used(selected))
    qaoa_energy = opt.fun + offset

    # Direct constraint check (also used above during post-selection, kept
    # here as a final confirmation on the chosen bitstring): not an
    # energy-tolerance proxy -- the tight penalty is an empirical,
    # QAOA-practical choice (validated on small test cases), NOT a rigorous
    # guarantee that it forbids every possible combination of simultaneous
    # violations on larger problems. Real hardware (run_ibm.py,
    # GCGCUUCGGCGC on ibm_marrakesh) demonstrated the plurality bitstring
    # CAN be infeasible even with this penalty -- that's why post-selection
    # exists above rather than trusting the raw top-shot bitstring blindly.
    is_feasible = not any(conflicting(q1, q2) for q1, q2 in combinations(selected, 2))

    return {
        "seq": seq, "n": len(seq), "num_qubits": result["num_qubits"],
        "qaoa_energy": qaoa_energy, "raw_stacking_energy": raw_energy,
        "structure": db, "shots_for_best_bitstring": best_shots,
        "total_shots": sum(counts.values()),
        "feasible": is_feasible, "penalty_used": penalty,
    }


if __name__ == "__main__":
    print("=== QAOA (Aer simulator, manual V2-primitive implementation) vs exact ===")
    for seq in ["GGGAAACCC", "GCGCUUCGGCGC"]:
        exact_energy, exact_selected = solve_exact(seq)
        exact_db = to_dot_bracket(len(seq), pairs_used(exact_selected))

        qaoa_out = run_qaoa_simulator(seq, reps=3, maxiter=300, restarts=8)

        print(f"seq={seq}  qubits={qaoa_out['num_qubits']}")
        print(f"  exact : energy={exact_energy:.2f}  structure={exact_db}")
        print(f"  QAOA  : energy={qaoa_out['qaoa_energy']:.2f}  "
              f"raw_stacking={qaoa_out['raw_stacking_energy']:.2f}  "
              f"structure={qaoa_out['structure']}  "
              f"best_bitstring_shots={qaoa_out['shots_for_best_bitstring']}/{qaoa_out['total_shots']}  "
              f"feasible={qaoa_out['feasible']}")
        print(f"  match: {'YES' if exact_db == qaoa_out['structure'] else 'NO'}")
        print()
