"""
BQM construction for the stacking-only RNA folding QUBO.

Variables: quartets (i,j) = stacked pair (i,j) on (i+1,j-1), real Turner2004
energy as linear bias (we MINIMIZE, so favorable stacks are negative bias).

Constraint encoding (penalty terms, quadratic):
  - Two quartets conflict (get a penalty edge) unless they either:
      (a) chain consistently -- q2's pair set and q1's pair set overlap
          exactly at one shared pair, forming one continuous stem, or
      (b) touch no common base index and don't cross.
  - Conflicting quartets get a positive quadratic penalty P, chosen larger
    than any achievable energy gain so the solver never pays it.

This is deliberately checked against the brute-force exact answer (from
validate_brute_force.py) on the same toy sequences before it's trusted on
anything larger or sent to hardware.
"""

import itertools
import dimod
from rna_qubo import build_quartets, quartets_crossing, quartet_bases, get_hairpin_energy


def implied_pairs(q):
    i, j = q
    return {(i, j), (i + 1, j - 1)}


def conflicting(q1, q2):
    if q1 == q2:
        return False
    p1, p2 = implied_pairs(q1), implied_pairs(q2)
    b1, b2 = quartet_bases(q1), quartet_bases(q2)
    shared_pairs = p1 & p2
    shared_bases = b1 & b2

    if shared_pairs:
        explained = set()
        for p in shared_pairs:
            explained |= {p[0], p[1]}
        if shared_bases <= explained:
            # clean chain (one stem extending by one stack) -- only a
            # problem if it somehow also crosses another selection
            return quartets_crossing(q1, q2)
        return True  # partial/inconsistent overlap -- not a valid chain

    if shared_bases:
        return True  # different stems fighting over the same base

    return quartets_crossing(q1, q2)


def build_bqm(seq, min_loop=3, penalty=None, include_hairpin=True):
    """
    include_hairpin=True (default) ports real hairpin-loop energy into the
    QUBO via a delta-correction technique:
      - every quartet q=(i,j) [inner pair p=(i+1,j-1)] gets a BASELINE
        linear addition of hp_energy(p), assuming q ends its stack chain
        here and p closes a hairpin.
      - if a continuation quartet next(q)=(i+1,j-1) also exists as a valid
        variable, a QUADRATIC correction term of -hp_energy(p) is added on
        (q, next(q)): cancels the wrongly-assumed hairpin baseline exactly
        when the chain actually continues past p.
    Requires NO auxiliary variables and NO new constraints, since quartets
    are already a linear chain (at most one possible continuation each) --
    unlike bulge/internal loops, which need auxiliary "at most one active
    child" constraints (real future work, not attempted here). Validated
    against dp_stack_hairpin.py, the exact ground truth for this same
    restricted model, before being trusted.
    include_hairpin=False reproduces the original stacking-only QUBO
    exactly, kept for comparison.
    """
    quartets = build_quartets(seq, min_loop)
    if not quartets:
        return dimod.BinaryQuadraticModel(vartype="BINARY"), quartets

    if penalty is None:
        # must exceed the maximum possible energy benefit from violating
        # a constraint, i.e. bigger than summing every favorable stack
        penalty = 2 * sum(abs(e) for e in quartets.values() if e < 0) + 10

    bqm = dimod.BinaryQuadraticModel(vartype="BINARY")
    for q, e in quartets.items():
        bqm.add_variable(q, e)

    if include_hairpin:
        for q in quartets:
            i, j = q
            p = (i + 1, j - 1)
            hp_e = get_hairpin_energy(seq, p[0], p[1])
            bqm.add_variable(q, hp_e)  # baseline: assume q ends the chain here
            next_q = (i + 1, j - 1)
            if next_q in quartets:
                bqm.add_interaction(q, next_q, -hp_e)  # cancel if chain continues

    for q1, q2 in itertools.combinations(quartets.keys(), 2):
        if conflicting(q1, q2):
            bqm.add_interaction(q1, q2, penalty)

    return bqm, quartets


def solve_exact(seq, min_loop=3):
    """dimod ExactSolver -- brute force over the BQM itself, used ONLY to
    cross-check against validate_brute_force.py's independent brute force.
    Same size limit applies (small sequences only)."""
    bqm, quartets = build_bqm(seq, min_loop)
    if len(bqm.variables) > 20:
        raise ValueError(f"{len(bqm.variables)} variables -- too many for ExactSolver")
    sampleset = dimod.ExactSolver().sample(bqm)
    best = sampleset.first
    selected = [q for q, v in best.sample.items() if v == 1]
    return best.energy, selected


if __name__ == "__main__":
    from validate_brute_force import brute_force_mfe, pairs_used, to_dot_bracket

    test_seqs = ["GGGAAACCC", "GCGCUUCGGCGC", "GGGGAAAACCCC", "GGAAUUCC", "CGCGCGAAAACGCGCG"]
    print(f"{'seq':<20}{'brute-force':<10}{'BQM exact':<10}{'match':<8}")
    for seq in test_seqs:
        bf_energy, bf_db, _ = brute_force_mfe(seq)
        bqm_energy, bqm_selected = solve_exact(seq)
        bqm_db = to_dot_bracket(len(seq), pairs_used(bqm_selected))
        match = "YES" if bf_db == bqm_db else "NO"
        print(f"{seq:<20}{bf_energy:<10.2f}{bqm_energy:<10.2f}{match:<8}")
        if match == "NO":
            print(f"    brute force: {bf_db}")
            print(f"    BQM        : {bqm_db}")
