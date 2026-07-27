"""
DP extending dp_stack_hairpin.py with real BULGE-loop energy (asymmetric
internal loop, gap on exactly one side, other side 0), bounded to a small
max bulge length. This is the ground truth the bulge QUBO port targets.

SCOPE, stated upfront: bulges only, not full two-sided internal loops.
A bulge has unpaired bases on exactly one side of the closing pair; a
full internal loop has unpaired bases on BOTH sides, which is a strictly
larger combinatorial space (paired lengths on each side, not just one)
left as documented future work. Multiloops remain unmodeled, same as
dp_stack_hairpin.py and dp_full_energy.py.

Max bulge length is bounded (default 4nt) for two reasons: real bulges
longer than a few nt are thermodynamically rare in favorable structures,
and bounding keeps both the DP and the eventual QUBO's variable count
tractable -- unbounded bulge length would need O(n) candidate successor
pairs per closing pair instead of O(max_bulge).
"""

import RNA
from rna_qubo import CANONICAL, get_stack_matrix, BP_TYPE


def dp_fold_bulge(seq, min_loop=3, max_bulge=4):
    n = len(seq)
    RNA.params_load_RNA_Turner2004()
    fc = RNA.fold_compound(seq)

    def hp(i0, j0):
        return fc.eval_hp_loop(i0 + 1, j0 + 1) / 100.0

    def il(i0, j0, k0, l0):
        return fc.eval_int_loop(i0 + 1, j0 + 1, k0 + 1, l0 + 1) / 100.0

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
                    cand = il(i, j, i + 1, j - 1) + V[i + 1][j - 1]
                    if cand < best_v:
                        best_v = cand
                        choice = ("bulge", i + 1, j - 1, 0, 0)

                for b in range(1, max_bulge + 1):
                    k, l = i + 1 + b, j - 1
                    if k >= l - min_loop:
                        break
                    if canon[k][l] and V[k][l] < NEG_INF:
                        cand = il(i, j, k, l) + V[k][l]
                        if cand < best_v:
                            best_v = cand
                            choice = ("bulge", k, l, b, 0)

                for b in range(1, max_bulge + 1):
                    k, l = i + 1, j - 1 - b
                    if k >= l - min_loop:
                        break
                    if canon[k][l] and V[k][l] < NEG_INF:
                        cand = il(i, j, k, l) + V[k][l]
                        if cand < best_v:
                            best_v = cand
                            choice = ("bulge", k, l, 0, b)

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
        choice = V_choice[i][j]
        if choice[0] == "bulge":
            _, k, l, b5, b3 = choice
            trace_v(k, l)

    trace_w(0, n - 1)
    return pairs


def to_dot_bracket(n, pairs):
    db = ["."] * n
    for (i, j) in pairs:
        db[i] = "("
        db[j] = ")"
    return "".join(db)


def bulge_mfe(seq, min_loop=3, max_bulge=4):
    W, V, V_choice, canon = dp_fold_bulge(seq, min_loop, max_bulge)
    pairs = traceback(seq, W, V, V_choice, canon, min_loop)
    return W[0][len(seq) - 1], to_dot_bracket(len(seq), pairs)


def compare_to_viennarna(seq):
    RNA.params_load_RNA_Turner2004()
    true_db, true_mfe = RNA.fold(seq)
    our_energy, our_db = bulge_mfe(seq)
    return {
        "seq": seq, "n": len(seq),
        "vienna_db": true_db, "vienna_mfe": true_mfe,
        "our_db": our_db, "our_energy": our_energy,
        "match": true_db == our_db,
    }


if __name__ == "__main__":
    print("=== self-check: recomputed energy vs ViennaRNA's own eval of OUR structure ===")
    test_seqs = ["GGGAAACCC", "GCGCUUCGGCGC", "GGGGAAAACCCC", "GGAAUUCC", "CGCGCGAAAACGCGCG"]
    for seq in test_seqs:
        our_e, our_db = bulge_mfe(seq)
        vienna_e_of_ours = RNA.energy_of_structure(seq, our_db, 0)
        vienna_db, vienna_mfe = RNA.fold(seq)
        agree = "OK" if abs(our_e - vienna_e_of_ours) < 0.05 else "MISMATCH"
        match = "MATCH" if our_db == vienna_db else "differ"
        print(f"{seq:<20} our={our_e:6.2f} vienna_eval={vienna_e_of_ours:6.2f} [{agree}]  "
              f"our_db={our_db}  vienna_mfe_db={vienna_db} [{match}]")
