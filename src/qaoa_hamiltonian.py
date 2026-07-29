"""
Same stacking-only QUBO as build_bqm.py, converted into a QAOA-ready Ising
Hamiltonian via qiskit-optimization. Kept as the SAME model, not a
reimplementation -- build_bqm() is imported directly, so both the D-Wave
and IBM paths run identical objective/constraint logic. Any difference in
results reflects the solver/hardware, not two different formulations.
"""

from qiskit_optimization import QuadraticProgram
from qiskit_optimization.translators import from_docplex_mp
from qiskit_optimization.converters import QuadraticProgramToQubo
from docplex.mp.model import Model

from build_bqm import build_bqm


def bqm_to_quadratic_program(bqm, quartets):
    """dimod BQM -> qiskit_optimization QuadraticProgram, same variables,
    same linear/quadratic coefficients, same penalty terms."""
    mdl = Model(name="rna_folding_stacking_only")
    var_list = list(bqm.variables)
    x = {v: mdl.binary_var(name=f"q_{v[0]}_{v[1]}") for v in var_list}

    linear = {x[v]: bqm.linear.get(v, 0.0) for v in var_list}
    obj = mdl.sum(coeff * x[v] for v, coeff in zip(var_list, [bqm.linear.get(v, 0.0) for v in var_list]))
    for (v1, v2), coeff in bqm.quadratic.items():
        obj += coeff * x[v1] * x[v2]

    mdl.minimize(obj)
    qp = from_docplex_mp(mdl)
    return qp, var_list, x


def build_qaoa_hamiltonian(seq, min_loop=3, penalty=None):
    bqm, quartets = build_bqm(seq, min_loop, penalty)
    qp, var_list, x = bqm_to_quadratic_program(bqm, quartets)

    conv = QuadraticProgramToQubo()
    qubo = conv.convert(qp)
    op, offset = qubo.to_ising()

    return {
        "quartets": quartets,
        "bqm": bqm,
        "var_list": var_list,
        "quadratic_program": qp,
        "qubo": qubo,
        "ising_op": op,
        "ising_offset": offset,
        "num_qubits": op.num_qubits,
    }


if __name__ == "__main__":
    for seq in ["GGGAAACCC", "GCGCUUCGGCGC"]:
        result = build_qaoa_hamiltonian(seq)
        print(f"{seq:<15} quartet variables: {len(result['var_list'])}   "
              f"qubits needed: {result['num_qubits']}   "
              f"Pauli terms: {len(result['ising_op'])}")
