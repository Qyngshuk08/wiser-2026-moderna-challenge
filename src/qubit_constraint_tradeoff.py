"""
Optional advanced task: "Analyze trade-offs between qubit count and
constraint enforcement."

Two things grow as sequence length (qubit count) increases:
  1. Constraint graph DENSITY -- already characterized in the D-Wave
     section (0.886 at n=20/15 vars down to 0.475 at n=100/540 vars).
  2. The GAP between two different penalty regimes:
       - "classical-exactness" penalty (build_bqm.py's default): sized to
         guarantee correctness against ARBITRARY constraint violations on
         classical exact solvers, scales with the SUM of every favorable
         energy in the problem.
       - "tight" penalty (qaoa_simulator.py's default): sized just above
         the SINGLE largest favorable energy, small enough for shallow
         QAOA to actually navigate the landscape.
  The classical penalty scales with problem size (sum over all quartets);
  the tight penalty does not (bounded by the single largest energy, which
  doesn't grow much with sequence length). The RATIO between them is a
  direct, quantifiable measure of how much harder constraint enforcement
  gets for NISQ-era solving as qubit count grows -- independent of and
  prior to any hardware noise effects.

Uses the SAME 6 sequences already used for D-Wave scaling
(scaling_results_run3_hybrid.json) so this is directly comparable to
existing findings, not a new arbitrary dataset.
"""

import json
from build_bqm import build_bqm

SEQUENCES = [
    ("GGACGGCGCUUCUACUCAAC", 20),
    ("GACGCUGGCAAAGAGCUCAUUUUGAACGAC", 30),
    ("CACUCGAGCUCUUUACGAAUUAAGCUUGCGGCACAGCUUA", 40),
    ("CAUCUCGAAUUAAAAUAGACAAAUUAGUAAGACAUACCGUAGAAGUCCGUUAUUCCAGUA", 60),
    ("GGCCGCUAGACACUUGCCCUCGGAAUCUUUGGGAGGAUGACUGGAAGCCGAGUAUAGGCAUAAAUAUCCUAGGAACGACU", 80),
    ("CAUCAAAUCCCAGGGCUUGCAUUCCCCACAACUUCCGUCAGUACUACGCCGAAGGCCAGAGGACCUUGAAAUUAGAAGCGGUGACUCUGUCAAAGUAGGC", 100),
]


def analyze(seq, n):
    # classical-exactness penalty, exactly as build_bqm.py computes it by
    # default (include_hairpin=True is the current default -- use it here
    # too, for consistency with the current state of the QUBO)
    bqm_default, quartets = build_bqm(seq, penalty=None)
    num_qubits = len(bqm_default.variables)
    num_edges = len(bqm_default.quadratic)
    max_possible_edges = num_qubits * (num_qubits - 1) // 2
    density = num_edges / max_possible_edges if max_possible_edges else 0.0

    classical_penalty = 2 * sum(abs(e) for e in quartets.values() if e < 0) + 10

    # tight penalty, exactly as qaoa_simulator.py computes it (calibrated
    # against the BQM's actual linear magnitudes, post-hairpin-port)
    preview_bqm, _ = build_bqm(seq, penalty=1.0)
    max_single = max((abs(v) for v in preview_bqm.linear.values()), default=1.0)
    tight_penalty = 1.5 * max_single + 1.0

    ratio = classical_penalty / tight_penalty if tight_penalty else float("inf")

    return {
        "seq": seq, "n": n, "num_qubits": num_qubits, "num_edges": num_edges,
        "density": round(density, 3),
        "classical_penalty": round(classical_penalty, 2),
        "tight_penalty": round(tight_penalty, 2),
        "ratio": round(ratio, 2),
    }


if __name__ == "__main__":
    results = [analyze(seq, n) for seq, n in SEQUENCES]

    print(f"{'n':>5}{'qubits':>8}{'density':>10}{'classical_penalty':>20}{'tight_penalty':>16}{'ratio':>9}")
    for r in results:
        print(f"{r['n']:>5}{r['num_qubits']:>8}{r['density']:>10.3f}"
              f"{r['classical_penalty']:>20.2f}{r['tight_penalty']:>16.2f}{r['ratio']:>9.2f}")

    with open("results/qubit_constraint_tradeoff_results.json", "w") as f:
        json.dump(results, f, indent=2)
