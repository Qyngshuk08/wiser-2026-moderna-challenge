"""
Exact brute-force validation: for small sequences, enumerate all feasible
subsets of quartets (respecting the one-pair-per-base and no-crossing
constraints), find the minimum-stacking-energy structure, and compare it
against ViennaRNA's true MFE structure.

This has to happen BEFORE any QAOA/annealing run. If the energy model
(stacking-only, real Turner04 values) doesn't get reasonably close to
ViennaRNA's MFE dot-bracket on toy sequences, no amount of quantum hardware
will fix that -- the model itself is wrong, not the solver.
"""

import itertools
import RNA
from rna_qubo import build_quartets, quartets_crossing, quartet_bases


def pairs_used(selected_quartets):
    """The set of actual base pairs (i,j) implied by a set of chained/selected
    quartets. Each quartet (i,j) implies pairs (i,j) and (i+1,j-1)."""
    pairs = set()
    for (i, j) in selected_quartets:
        pairs.add((i, j))
        pairs.add((i + 1, j - 1))
    return pairs


def feasible(selected_quartets):
    """Check one-pair-per-base and no-pseudoknot constraints."""
    pairs = pairs_used(selected_quartets)
    seen_bases = {}
    for (i, j) in pairs:
        for b in (i, j):
            if b in seen_bases and seen_bases[b] != (i, j):
                return False
            seen_bases[b] = (i, j)
    for q1, q2 in itertools.combinations(selected_quartets, 2):
        if quartets_crossing(q1, q2):
            return False
    return True


def to_dot_bracket(n, pairs):
    db = ["."] * n
    for (i, j) in pairs:
        db[i] = "("
        db[j] = ")"
    return "".join(db)


def brute_force_mfe(seq, min_loop=3):
    quartets = build_quartets(seq, min_loop)
    keys = list(quartets.keys())
    if len(keys) > 20:
        raise ValueError(f"{len(keys)} quartets -- too many for brute force, use a solver instead")

    best_energy = 0.0  # empty structure = 0 energy, always feasible baseline
    best_selection = ()
    for r in range(1, len(keys) + 1):
        for combo in itertools.combinations(keys, r):
            if feasible(combo):
                e = sum(quartets[q] for q in combo)
                if e < best_energy:
                    best_energy = e
                    best_selection = combo
    pairs = pairs_used(best_selection)
    return best_energy, to_dot_bracket(len(seq), pairs), best_selection


def compare_to_viennarna(seq):
    RNA.params_load_RNA_Turner2004()
    true_db, true_mfe = RNA.fold(seq)
    our_energy, our_db, _ = brute_force_mfe(seq)
    print(f"seq               : {seq}")
    print(f"ViennaRNA MFE     : {true_db}  ({true_mfe:.2f} kcal/mol, full model)")
    print(f"stacking-only QUBO: {our_db}  ({our_energy:.2f} kcal/mol, stacking-only)")
    print(f"match: {'YES' if true_db == our_db else 'NO'}")
    print()


if __name__ == "__main__":
    test_seqs = [
        "GGGAAACCC",        # simple hairpin, should be trivial
        "GCGCUUCGGCGC",     # slightly longer hairpin-forming sequence
        "GGGGAAAACCCC",     # longer stem
        "GGAAUUCC",         # short, weak stem
        "CGCGCGAAAACGCGCG", # two possible stems, tests optimality
    ]
    for s in test_seqs:
        compare_to_viennarna(s)
