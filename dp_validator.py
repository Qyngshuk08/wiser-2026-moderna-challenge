"""
Exact O(n^3) DP for the stacking-only energy model, replacing
validate_brute_force.py's ~20-quartet ceiling.

This is a Zuker/Nussinov-style two-table DP (W = best over a window,
V = best given the two endpoints are paired), specialized to a model with
ONLY stacking energy and no loop penalties (since that's the model we're
actually solving with the QUBO). Because there's no hairpin/bulge/internal-
loop penalty, closing a pair with nothing enclosed costs 0 -- so V(i,j)
either continues a stack inward (real ΔG) or terminates for free.

This DP computes the EXACT global optimum of the exact same stacking-only
model the QUBO approximates, at O(n^3) instead of exponential brute force.
It is NOT computing ViennaRNA's real MFE (that already exists via RNA.fold);
it's giving us a fast, exact ground truth for OUR simplified model, so we
can check the model's match rate against real MFE at real scale instead of
just five toy sequences.
"""

import RNA
from rna_qubo import CANONICAL, BP_TYPE, get_stack_matrix


def dp_fold(seq, min_loop=3):
    n = len(seq)
    stack = get_stack_matrix()
    NEG_INF = float("inf")

    canon = [[False] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if (seq[i], seq[j]) in CANONICAL:
                canon[i][j] = True

    def stack_energy(i, j):
        # stacking energy of pair (i,j) on inner pair (i+1,j-1)
        bt_o = BP_TYPE[(seq[i], seq[j])]
        bt_i = BP_TYPE[(seq[i + 1], seq[j - 1])]
        return stack[bt_o][bt_i] / 100.0

    # W[i][j]: best (most negative) energy for window [i, j], i,j not
    # necessarily paired. V[i][j]: best energy given (i,j) ARE paired
    # (only defined where canon[i][j] holds).
    W = [[0.0] * n for _ in range(n)]
    V = [[NEG_INF] * n for _ in range(n)]

    for span in range(min_loop + 1, n):
        for i in range(0, n - span):
            j = i + span
            # V(i,j): only if (i,j) can pair
            if canon[i][j]:
                best_v = 0.0  # terminate here, no loop penalty in this model
                if canon[i + 1][j - 1] and (j - 1) - (i + 1) >= min_loop + 1 - 2:
                    if j - 1 > i + 1:  # inner window must be a valid pair site
                        cand = stack_energy(i, j) + V[i + 1][j - 1]
                        if cand < best_v:
                            best_v = cand
                V[i][j] = best_v

            best_w = W[i + 1][j]  # i unpaired
            if W[i][j - 1] < best_w:
                best_w = W[i][j - 1]  # j unpaired
            for k in range(i, j - min_loop):
                if canon[k][j] and V[k][j] < NEG_INF:
                    left = W[i][k - 1] if k > i else 0.0
                    cand = left + V[k][j]
                    if cand < best_w:
                        best_w = cand
            W[i][j] = best_w

    return W, V, canon


def traceback(seq, W, V, canon, min_loop=3):
    n = len(seq)
    pairs = []

    def trace_w(i, j):
        if j - i <= min_loop:
            return
        if W[i][j] == W[i + 1][j]:
            trace_w(i + 1, j)
            return
        if W[i][j] == W[i][j - 1]:
            trace_w(i, j - 1)
            return
        for k in range(i, j - min_loop):
            if canon[k][j]:
                left = W[i][k - 1] if k > i else 0.0
                if abs(W[i][j] - (left + V[k][j])) < 1e-9:
                    if k > i:
                        trace_w(i, k - 1)
                    trace_v(k, j)
                    return

    def trace_v(i, j):
        pairs.append((i, j))
        if canon[i + 1][j - 1] if (j - 1 > i + 1 or j - 1 >= i + 1) else False:
            pass
        if j - 1 > i + 1 and canon[i + 1][j - 1]:
            from rna_qubo import BP_TYPE
            stack = get_stack_matrix()
            bt_o = BP_TYPE[(seq[i], seq[j])]
            bt_i = BP_TYPE[(seq[i + 1], seq[j - 1])]
            e = stack[bt_o][bt_i] / 100.0
            if abs(V[i][j] - (e + V[i + 1][j - 1])) < 1e-9:
                trace_v(i + 1, j - 1)

    trace_w(0, n - 1)
    return pairs


def to_dot_bracket(n, pairs):
    db = ["."] * n
    for (i, j) in pairs:
        db[i] = "("
        db[j] = ")"
    return "".join(db)


def dp_mfe(seq, min_loop=3):
    W, V, canon = dp_fold(seq, min_loop)
    pairs = traceback(seq, W, V, canon, min_loop)
    energy = W[0][len(seq) - 1]
    return energy, to_dot_bracket(len(seq), pairs)


def compare_to_viennarna(seq):
    RNA.params_load_RNA_Turner2004()
    true_db, true_mfe = RNA.fold(seq)
    our_energy, our_db = dp_mfe(seq)
    match = true_db == our_db
    return {
        "seq": seq, "n": len(seq),
        "vienna_db": true_db, "vienna_mfe": true_mfe,
        "our_db": our_db, "our_energy": our_energy,
        "match": match,
    }


if __name__ == "__main__":
    from validate_brute_force import compare_to_viennarna as brute_compare, brute_force_mfe

    # cross-check against the independent brute-force solver on the same
    # toy sequences before trusting the DP on anything new
    print("=== cross-check DP vs brute force (must agree exactly) ===")
    for seq in ["GGGAAACCC", "GCGCUUCGGCGC", "GGGGAAAACCCC", "GGAAUUCC", "CGCGCGAAAACGCGCG"]:
        bf_e, bf_db, _ = brute_force_mfe(seq)
        dp_e, dp_db = dp_mfe(seq)
        ok = "OK" if (abs(bf_e - dp_e) < 1e-6 and bf_db == dp_db) else "MISMATCH"
        print(f"{seq:<20} brute={bf_e:.2f}/{bf_db}   dp={dp_e:.2f}/{dp_db}   {ok}")
