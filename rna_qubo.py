"""
QUBO formulation for RNA secondary structure (MFE) prediction using REAL
Turner2004 nearest-neighbor stacking energies pulled directly from ViennaRNA,
rather than the tunable heuristic constants used in prior quantum-annealing
RNA folding papers (Fox et al. 2022, Zaborniak et al. 2022).

Variables: quartets q(i,j) representing a stacked pair -- base pair (i,j)
immediately stacked on base pair (i+1,j-1). This mirrors the MIS-style
formulation used in recent gate-based work (arXiv:2505.05782) but replaces
their tunable reward constants with actual ΔG values.

Limitations (documented, not hidden):
  - Only captures stacking energy. Hairpin/bulge/internal-loop/multiloop
    entropic terms are NOT included in this first version (they are not
    naturally quadratic -- length-dependent). This means MFE will be
    approximated by "maximize favorable stacking," not the true Turner04
    free energy. Loop terms are the planned second-pass extension.
  - No pseudoknots (crossing pairs forbidden by construction).
"""

import itertools
import RNA

CANONICAL = {("A", "U"), ("U", "A"), ("C", "G"), ("G", "C"), ("G", "U"), ("U", "G")}

# ViennaRNA base-pair type indices (bptype), used to index the stack[][] matrix.
# 0 = none/invalid, 1=CG, 2=GC, 3=GU, 4=UG, 5=AU, 6=UA, 7=non-standard
BP_TYPE = {
    ("C", "G"): 1, ("G", "C"): 2,
    ("G", "U"): 3, ("U", "G"): 4,
    ("A", "U"): 5, ("U", "A"): 6,
}


def get_stack_matrix():
    """Real Turner2004 stacking energies (10 cal/mol units) from ViennaRNA."""
    RNA.params_load_RNA_Turner2004()
    p = RNA.param()
    return p.stack  # 8x8 int matrix


def valid_pairs(seq, min_loop=3):
    """All (i, j), 0-indexed, i < j, canonical Watson-Crick/wobble, respecting
    the minimum hairpin loop length (ViennaRNA default = 3)."""
    n = len(seq)
    pairs = []
    for i in range(n):
        for j in range(i + min_loop + 1, n):
            if (seq[i], seq[j]) in CANONICAL:
                pairs.append((i, j))
    return pairs


def build_quartets(seq, min_loop=3):
    """Enumerate stacked-pair variables and their real ΔG stacking energy.

    A quartet is valid if (i,j) and (i+1,j-1) are both canonical pairs.
    Energy convention: stack[][] values are already negative (favorable)
    for stable stacks, in units of 10 cal/mol. We convert to kcal/mol.
    """
    stack = get_stack_matrix()
    pset = set(valid_pairs(seq, min_loop))
    quartets = {}  # key: (i,j) meaning pair(i,j) stacked on pair(i+1,j-1)
    for (i, j) in pset:
        inner = (i + 1, j - 1)
        if inner in pset:
            bt_outer = BP_TYPE[(seq[i], seq[j])]
            bt_inner = BP_TYPE[(seq[i + 1], seq[j - 1])]
            e = stack[bt_outer][bt_inner] / 100.0  # 10cal/mol -> kcal/mol
            quartets[(i, j)] = e
    return quartets


def quartet_bases(q):
    """The 4 sequence positions a quartet touches: i, j, i+1, j-1."""
    i, j = q
    return {i, j, i + 1, j - 1}


def quartets_crossing(q1, q2):
    """True if the two stacked pairs would form a pseudoknot."""
    i1, j1 = q1
    i2, j2 = q2
    # crossing if exactly one of i2,j2 lies strictly inside (i1,j1)
    def crosses(a, b, c, d):
        return (a < c < b < d) or (c < a < d < b)
    return crosses(i1, j1, i2, j2)


def get_hairpin_energy(seq, i, j):
    """Real hairpin-loop closing energy for pair (i,j), 0-indexed, via
    ViennaRNA's own evaluator (not hand-indexed tables) -- same validated
    approach as dp_full_energy.py / dp_stack_hairpin.py. i,j must be a
    valid canonical pair with a legal loop length."""
    fc = get_hairpin_energy._fc_cache.get(seq)
    if fc is None:
        RNA.params_load_RNA_Turner2004()
        fc = RNA.fold_compound(seq)
        get_hairpin_energy._fc_cache = {seq: fc}  # single-entry cache, one seq at a time
    return fc.eval_hp_loop(i + 1, j + 1) / 100.0


get_hairpin_energy._fc_cache = {}


if __name__ == "__main__":
    seq = "GGGAAACCC"
    q = build_quartets(seq)
    print(f"seq={seq}  n={len(seq)}  #quartets={len(q)}")
    for k, v in sorted(q.items()):
        print(f"  stack at {k} (outer pair {k}, inner pair ({k[0]+1},{k[1]-1})): {v:.2f} kcal/mol")
