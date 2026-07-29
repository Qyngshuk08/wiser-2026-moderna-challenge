"""
DP for the STACKING + HAIRPIN model only (no bulge/internal loops, no
multiloops). This is deliberately the exact energy model the QUBO
extension in rna_qubo.py/build_bqm.py targets -- quartets already only
support strict adjacent stacking (i,j)-(i+1,j-1), so bulge/internal loops
(which require a flexible gap) are out of scope for that formulation
without auxiliary "single active child" constraint variables (a harder,
separate task -- see README).

Built first, before writing any QUBO code, to answer honestly: how much of
dp_full_energy.py's 10.0% -> 53.4% improvement is attributable to hairpin
energy alone (which the QUBO CAN capture) vs bulge/internal loops (which it
currently CANNOT)? That answer determines whether the QUBO port is worth
doing on its own or needs the harder bulge/internal work to matter.
"""

import RNA
from rna_qubo import CANONICAL, get_stack_matrix, BP_TYPE


def dp_fold_stack_hairpin(seq, min_loop=3):
    n = len(seq)
    stack = get_stack_matrix()
    RNA.params_load_RNA_Turner2004()
    fc = RNA.fold_compound(seq)

    def hp(i0, j0):
        return fc.eval_hp_loop(i0 + 1, j0 + 1) / 100.0

    def stack_e(i0, j0):
        bt_o = BP_TYPE[(seq[i0], seq[j0])]
        bt_i = BP_TYPE[(seq[j0 - 1], seq[i0 + 1])]  # REVERSED -- see rna_qubo.py fix note
        return stack[bt_o][bt_i] / 100.0

    canon = [[False] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if (seq[i], seq[j]) in CANONICAL:
                canon[i][j] = True

    NEG_INF = float("inf")
    W = [[0.0] * n for _ in range(n)]
    V = [[NEG_INF] * n for _ in range(n)]
    V_choice = [[None] * n for _ in range(n)]

    for span in range(min_loop + 1, n):
        for i in range(0, n - span):
            j = i + span
            if canon[i][j]:
                best_v = hp(i, j)
                choice = ("hairpin",)
                if j - 1 > i + 1 and canon[i + 1][j - 1] and V[i + 1][j - 1] < NEG_INF:
                    cand = stack_e(i, j) + V[i + 1][j - 1]
                    if cand < best_v:
                        best_v = cand
                        choice = ("stack",)
                V[i][j] = best_v
                V_choice[i][j] = choice

            best_w = W[i + 1][j]
            if W[i][j - 1] < best_w:
                best_w = W[i][j - 1]
            for k in range(i, j - min_loop):
                if canon[k][j] and V[k][j] < NEG_INF:
                    left = W[i][k - 1] if k > i else 0.0
                    cand = left + V[k][j]
                    if cand < best_w:
                        best_w = cand
            W[i][j] = best_w

    return W, V, V_choice, canon


def traceback(seq, W, V, V_choice, canon, min_loop=3):
    n = len(seq)
    pairs = []

    def trace_w(i, j):
        if j - i <= min_loop:
            return
        if abs(W[i][j] - W[i + 1][j]) < 1e-9:
            trace_w(i + 1, j)
            return
        if abs(W[i][j] - W[i][j - 1]) < 1e-9:
            trace_w(i, j - 1)
            return
        for k in range(i, j - min_loop):
            if canon[k][j]:
                left = W[i][k - 1] if k > i else 0.0
                if abs(W[i][j] - (left + V[k][j])) < 1e-7:
                    if k > i:
                        trace_w(i, k - 1)
                    trace_v(k, j)
                    return

    def trace_v(i, j):
        pairs.append((i, j))
        if V_choice[i][j][0] == "stack":
            trace_v(i + 1, j - 1)

    trace_w(0, n - 1)
    return pairs


def to_dot_bracket(n, pairs):
    db = ["."] * n
    for (i, j) in pairs:
        db[i] = "("
        db[j] = ")"
    return "".join(db)


def stack_hairpin_mfe(seq, min_loop=3):
    W, V, V_choice, canon = dp_fold_stack_hairpin(seq, min_loop)
    pairs = traceback(seq, W, V, V_choice, canon, min_loop)
    return W[0][len(seq) - 1], to_dot_bracket(len(seq), pairs)


def compare_to_viennarna(seq):
    RNA.params_load_RNA_Turner2004()
    true_db, true_mfe = RNA.fold(seq)
    our_energy, our_db = stack_hairpin_mfe(seq)
    return {
        "seq": seq, "n": len(seq),
        "vienna_db": true_db, "vienna_mfe": true_mfe,
        "our_db": our_db, "our_energy": our_energy,
        "match": true_db == our_db,
    }


if __name__ == "__main__":
    print("=== self-check: our energy vs ViennaRNA's own eval of OUR structure ===")
    for seq in ["GGGAAACCC", "GCGCUUCGGCGC", "GGGGAAAACCCC", "GGAAUUCC", "CGCGCGAAAACGCGCG"]:
        e, db = stack_hairpin_mfe(seq)
        vienna_e = RNA.energy_of_structure(seq, db, 0)
        vienna_db, vienna_mfe = RNA.fold(seq)
        agree = "OK" if abs(e - vienna_e) < 0.05 else "MISMATCH"
        match = "MATCH" if db == vienna_db else "differ"
        print(f"{seq:<20} our={e:6.2f} vienna_eval={vienna_e:6.2f} [{agree}]  our_db={db}  vienna_db={vienna_db} [{match}]")
