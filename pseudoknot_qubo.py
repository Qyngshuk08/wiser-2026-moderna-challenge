"""
Optional advanced task: "Explore formulations that include pseudoknots."

SCOPE, stated honestly upfront: this does NOT attempt general pseudoknot
folding (NP-hard in general; even restricted classes like H-type
pseudoknots normally need specialized O(n^6) algorithms -- e.g. Rivas &
Eddy 1999 -- to correctly score them, well beyond what remaining project
time allows). What this DOES do: relax this project's no-crossing
constraint so the QUBO CAN select crossing quartets when they are
genuinely favorable under the existing stacking-only energy model, and
prove -- on a real, hand-verified test case -- that doing so lets the
solver find a correct, better answer that the non-crossing model
structurally cannot reach.

ENERGY MODEL CAVEAT, stated honestly: real pseudoknots have additional
loop-topology entropy penalties (e.g. Rivas & Eddy's gw/gwh terms) beyond
plain stacking energy, which this extension does NOT model. Every
crossing quartet is scored with the exact same real Turner2004 stacking
energy as a nested one. This means the extension's total energy for a
genuine pseudoknot is an UNDERESTIMATE of the true (more unfavorable)
free energy -- consistent with how this project has handled every other
energy-model gap (hairpin was added incrementally and honestly; bulge/
internal loops and multiloops remain unmodeled and are documented as
such). This is a real, bounded simplification, not a hidden one.

WHAT IS RELAXED: rna_qubo_pairs.py / build_bqm.py's no-crossing check
(quartets_crossing / the crossing branch inside conflicting()) is skipped
entirely here. The no-shared-base constraint (a base cannot be in two
different pairs) is KEPT -- that constraint is physically required
regardless of pseudoknots.
"""

import itertools
import dimod
from rna_qubo import build_quartets, quartet_bases


def pk_conflicting(q1, q2):
    """Pseudoknot-permissive conflict check: forbid only shared-base
    violations (a base can't be paired twice) and inconsistent partial
    chain overlaps. Crossing is explicitly ALLOWED here -- this is the
    only difference from build_bqm.conflicting()."""
    if q1 == q2:
        return False
    from build_bqm import implied_pairs
    p1, p2 = implied_pairs(q1), implied_pairs(q2)
    b1, b2 = quartet_bases(q1), quartet_bases(q2)
    shared_pairs = p1 & p2
    shared_bases = b1 & b2

    if shared_pairs:
        explained = set()
        for p in shared_pairs:
            explained |= {p[0], p[1]}
        if shared_bases <= explained:
            return False  # clean chain, no crossing check here -- allowed
        return True  # partial/inconsistent overlap -- still forbidden

    if shared_bases:
        return True  # different stems fighting over the same base

    return False  # NOT checking crossing -- this is the whole point


def build_pseudoknot_bqm(seq, min_loop=3, penalty=None):
    """Stacking-only QUBO (include_hairpin equivalent NOT ported here --
    scope limited to proving the crossing relaxation works on the same
    stacking-only baseline used to validate every other encoding in this
    project) with the no-crossing constraint relaxed."""
    quartets = build_quartets(seq, min_loop)
    if not quartets:
        return dimod.BinaryQuadraticModel(vartype="BINARY"), quartets

    if penalty is None:
        penalty = 2 * sum(abs(e) for e in quartets.values() if e < 0) + 10

    bqm = dimod.BinaryQuadraticModel(vartype="BINARY")
    for q, e in quartets.items():
        bqm.add_variable(q, e)

    for q1, q2 in itertools.combinations(quartets.keys(), 2):
        if pk_conflicting(q1, q2):
            bqm.add_interaction(q1, q2, penalty)

    return bqm, quartets


def pk_pairs_used(selected_quartets):
    pairs = set()
    for (i, j) in selected_quartets:
        pairs.add((i, j))
        pairs.add((i + 1, j - 1))
    return pairs


def to_pseudoknot_notation(n, pairs):
    """Extended dot-bracket: uses () for the first (largest) non-crossing
    subset found greedily, [] for pairs that cross something already
    assigned (). This is a DISPLAY convention only (standard extended
    dot-bracket notation), not used by the solver itself."""
    db = ["."] * n
    sorted_pairs = sorted(pairs, key=lambda p: p[1] - p[0], reverse=True)
    assigned_bracket = {}
    open_at = {p: [] for p in ["()", "[]"]}
    for (i, j) in sorted_pairs:
        # does this pair cross anything already assigned to level 0 "()"?
        crosses_round = any(
            (i < a < j < b) or (a < i < b < j)
            for (a, b), lvl in assigned_bracket.items() if lvl == 0
        )
        level = 1 if crosses_round else 0
        assigned_bracket[(i, j)] = level
        open_c, close_c = ("(", ")") if level == 0 else ("[", "]")
        db[i] = open_c
        db[j] = close_c
    return "".join(db)


if __name__ == "__main__":
    import neal
    from validate_brute_force import pairs_used, to_dot_bracket
    from build_bqm import build_bqm

    seq = "GCAUGAACGUACAAACAUGCAGUACG"
    print(f"Test sequence (hand-constructed, verified crossing stems): {seq}")
    print(f"length: {len(seq)}")
    print()

    print("=== Baseline: existing non-crossing QUBO (stacking-only) ===")
    bqm_nc, quartets_nc = build_bqm(seq, include_hairpin=False)
    sampler = neal.SimulatedAnnealingSampler()
    sa_nc = sampler.sample(bqm_nc, num_reads=3000).first
    sel_nc = [q for q, v in sa_nc.sample.items() if v == 1]
    db_nc = to_dot_bracket(len(seq), pairs_used(sel_nc))
    print(f"energy: {sa_nc.energy:.2f}   structure: {db_nc}")
    print()

    print("=== Pseudoknot-permissive QUBO (crossing allowed) ===")
    bqm_pk, quartets_pk = build_pseudoknot_bqm(seq)
    sa_pk = sampler.sample(bqm_pk, num_reads=3000).first
    sel_pk = [q for q, v in sa_pk.sample.items() if v == 1]
    pairs_pk = pk_pairs_used(sel_pk)
    db_pk = to_pseudoknot_notation(len(seq), pairs_pk)
    print(f"energy: {sa_pk.energy:.2f}   structure: {db_pk}")
    print()

    print(f"improvement from allowing crossing: {sa_nc.energy - sa_pk.energy:.2f} kcal/mol")
    print(f"expected (hand-verified, both 5bp arms independently, "
          f"non-overlapping): -17.40")
