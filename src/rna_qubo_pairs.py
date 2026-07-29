"""
Optional advanced task: "Compare multiple quantum encodings."

ENCODING A (this project's primary formulation, rna_qubo.py/build_bqm.py):
  Variables = QUARTETS -- one binary variable per stacked pair (i,j) on
  (i+1,j-1). Stacking energy is baked directly into the variable's linear
  bias. A "chain" of quartets naturally represents a stem; the "one pair
  per base" and "no crossing" constraints are penalty edges between
  quartets.

ENCODING B (this file): the RAW PAIR encoding used in prior literature
  (Fox, DePrince, Skolnick 2022; Zaborniak et al. 2022) -- one binary
  variable per valid base pair (i,j), full stop. A pair alone gets ZERO
  linear reward (a lone pair has no stacking, so no reward under a
  stacking-only energy model); stacking energy is a QUADRATIC bonus
  between two variables representing adjacent pairs (i,j) and (i+1,j-1).
  Constraints (one pair per base, no crossing) are penalty edges between
  ALL conflicting pairs, same idea as Encoding A but over a different,
  larger variable set.

These encode the IDENTICAL physical model (stacking-only MFE) two
different ways. Comparing them means comparing: qubit count, constraint
graph density, and whether both actually reach the same ground truth --
not comparing two different physics models.
"""

import itertools
import dimod
from rna_qubo import valid_pairs, get_stack_matrix, BP_TYPE


def build_pair_bqm(seq, min_loop=3, penalty=None):
    stack = get_stack_matrix()
    pairs = valid_pairs(seq, min_loop)
    pair_set = set(pairs)

    if not pairs:
        return dimod.BinaryQuadraticModel(vartype="BINARY"), {}

    bqm = dimod.BinaryQuadraticModel(vartype="BINARY")
    for p in pairs:
        bqm.add_variable(p, 0.0)  # a lone pair gets zero reward, same as the real model

    stack_terms = {}
    for (i, j) in pairs:
        inner = (i + 1, j - 1)
        if inner in pair_set:
            bt_o = BP_TYPE[(seq[i], seq[j])]
            bt_i = BP_TYPE[(seq[j - 1], seq[i + 1])]  # REVERSED -- see rna_qubo.py fix note
            e = stack[bt_o][bt_i] / 100.0
            bqm.add_interaction((i, j), inner, e)
            stack_terms[((i, j), inner)] = e

    if penalty is None:
        penalty = 2 * sum(abs(e) for e in stack_terms.values() if e < 0) + 10

    def conflicting_pairs(p1, p2):
        i1, j1 = p1
        i2, j2 = p2
        if p1 == p2:
            return False
        shared_base = len({i1, j1} & {i2, j2}) > 0
        if shared_base:
            return True
        # crossing check
        a, b, c, d = i1, j1, i2, j2
        return (a < c < b < d) or (c < a < d < b)

    for p1, p2 in itertools.combinations(pairs, 2):
        if conflicting_pairs(p1, p2):
            bqm.add_interaction(p1, p2, penalty)

    return bqm, stack_terms


def solve_exact_pairs(seq, min_loop=3):
    bqm, stack_terms = build_pair_bqm(seq, min_loop)
    if len(bqm.variables) > 20:
        raise ValueError(f"{len(bqm.variables)} variables -- too many for ExactSolver")
    sampleset = dimod.ExactSolver().sample(bqm)
    best = sampleset.first
    selected = [p for p, v in best.sample.items() if v == 1]
    return best.energy, selected


def pairs_to_dot_bracket(n, selected_pairs):
    db = ["."] * n
    for (i, j) in selected_pairs:
        db[i] = "("
        db[j] = ")"
    return "".join(db)


if __name__ == "__main__":
    from build_bqm import solve_exact as solve_exact_quartets
    from validate_brute_force import pairs_used, to_dot_bracket

    test_seqs = ["GGGAAACCC", "GCGCUUCGGCGC", "GGGGAAAACCCC", "GGAAUUCC", "CGCGCGAAAACGCGCG"]
    print("Both encodings compared on the SAME stacking-only model "
          "(include_hairpin=False on the quartet side) for a fair,\n"
          "like-for-like comparison -- Encoding B has no hairpin term.\n")
    print(f"{'seq':<20}{'A:quartets qubits':>18}{'B:pairs qubits':>16}{'A energy':>11}{'B energy':>11}{'match':>8}")
    for seq in test_seqs:
        e_a, sel_a = solve_exact_quartets(seq, include_hairpin=False)
        db_a = to_dot_bracket(len(seq), pairs_used(sel_a))

        e_b, sel_b = solve_exact_pairs(seq)
        db_b = pairs_to_dot_bracket(len(seq), sel_b)

        from build_bqm import build_bqm
        bqm_a, _ = build_bqm(seq, include_hairpin=False)
        bqm_b, _ = build_pair_bqm(seq)

        match = "YES" if db_a == db_b else "NO"
        print(f"{seq:<20}{len(bqm_a.variables):>18}{len(bqm_b.variables):>16}{e_a:>11.2f}{e_b:>11.2f}{match:>8}")
        if match == "NO":
            print(f"    A (quartets): {db_a}")
            print(f"    B (pairs):    {db_b}")
