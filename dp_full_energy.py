"""
DP with REAL hairpin, bulge, and internal loop energies, pulled directly
from ViennaRNA's own evaluator functions (fold_compound.eval_hp_loop,
fold_compound.eval_int_loop) rather than hand-indexed parameter tables.
This eliminates indexing/unit mistakes by construction: eval_int_loop with
adjacent inner pair (k=i+1, l=j-1) IS the real stacking energy, verified
above to reconstruct ViennaRNA's true total energy on a known structure
exactly (-1.20 kcal/mol on GGGAAACCC, bit for bit).

What this DOES model, with real Turner2004 numbers:
  - stacking (as a special case of internal loop with 0x0 unpaired)
  - bulge loops (internal loop, one side 0 unpaired)
  - internal loops (both sides >=1 unpaired), including asymmetry and
    terminal mismatch energies -- all computed by ViennaRNA itself, not
    reimplemented

What this does NOT model (documented, not hidden):
  - Multiloops. Any pair enclosing more than one independent branch still
    falls back to a zero-cost closing (W(i+1,j-1)), same simplification as
    rna_qubo.py / dp_validator.py. This is a real, acknowledged gap -- a
    correct multiloop model (MLbase/MLintern/MLclosing with a WM table)
    is future work, out of scope for this validation pass.

Complexity: O(n^2) states, O(n^2) work per V(i,j) state (scanning inner
pair candidates) = O(n^4) total, further capped by max_loop_size (mirrors
ViennaRNA's own default internal-loop-size cutoff of 30 for tractability).
This is a VALIDATOR, not the production QUBO solver -- it does not need to
scale past ~60nt to be useful; it exists to check whether a full-energy
QUBO extension would be worth building before committing days to it.
"""

import RNA
from rna_qubo import CANONICAL


def dp_fold_full(seq, min_loop=3, max_loop_size=30):
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
    V_choice = [[None] * n for _ in range(n)]  # for traceback

    for span in range(min_loop + 1, n):
        for i in range(0, n - span):
            j = i + span
            if canon[i][j]:
                best_v = hp(i, j)
                choice = ("hairpin",)
                # scan inner pairs (k,l), i<k<l<j, bounded loop size
                for k in range(i + 1, j):
                    max_l = min(j - 1, k + (max_loop_size - (k - i - 1)) + 1)
                    for l in range(k + 1, min(j, max_l + 1)):
                        unpaired = (k - i - 1) + (j - l - 1)
                        if unpaired > max_loop_size:
                            continue
                        if canon[k][l] and V[k][l] < NEG_INF:
                            cand = il(i, j, k, l) + V[k][l]
                            if cand < best_v:
                                best_v = cand
                                choice = ("loop", k, l)
                # NOTE: no multiloop fallback. An earlier version allowed a
                # free zero-cost "close with nothing enclosed" option here,
                # intended to represent multi-branch loops. Since 0 is always
                # less than any real hairpin penalty, that path was being
                # selected on EVERY hairpin, not just genuine multiloops --
                # fabricating free structure rather than approximating one.
                # Removed. This DP currently supports NO multiloops at all
                # (a conservative gap, not a silent cheat): structures that
                # require a real 3+-way branch will not be found. A correct
                # fix needs a separate WM table with MLbase/MLintern/
                # MLclosing costs and branch-count tracking -- out of scope
                # for this pass, tracked as a real limitation below.
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


def traceback_full(seq, W, V, V_choice, canon, min_loop=3):
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
        if choice[0] == "loop":
            _, k, l = choice
            trace_v(k, l)
        # "hairpin": nothing further inside

    trace_w(0, n - 1)
    return pairs


def to_dot_bracket(n, pairs):
    db = ["."] * n
    for (i, j) in pairs:
        db[i] = "("
        db[j] = ")"
    return "".join(db)


def full_mfe(seq, min_loop=3, max_loop_size=30):
    W, V, V_choice, canon = dp_fold_full(seq, min_loop, max_loop_size)
    pairs = traceback_full(seq, W, V, V_choice, canon, min_loop)
    return W[0][len(seq) - 1], to_dot_bracket(len(seq), pairs)


def compare_to_viennarna(seq):
    RNA.params_load_RNA_Turner2004()
    true_db, true_mfe = RNA.fold(seq)
    our_energy, our_db = full_mfe(seq)
    return {
        "seq": seq, "n": len(seq),
        "vienna_db": true_db, "vienna_mfe": true_mfe,
        "our_db": our_db, "our_energy": our_energy,
        "match": true_db == our_db,
    }


if __name__ == "__main__":
    print("=== sanity check: does our recomputed energy match ViennaRNA's own "
          "eval_structure on OUR predicted structure? (validates the energy "
          "function itself, independent of whether the fold matches) ===")
    test_seqs = ["GGGAAACCC", "GCGCUUCGGCGC", "GGGGAAAACCCC", "GGAAUUCC", "CGCGCGAAAACGCGCG"]
    for seq in test_seqs:
        RNA.params_load_RNA_Turner2004()
        our_e, our_db = full_mfe(seq)
        vienna_e_of_our_structure = RNA.energy_of_structure(seq, our_db, 0)
        vienna_db, vienna_mfe = RNA.fold(seq)
        agree = "OK" if abs(our_e - vienna_e_of_our_structure) < 0.05 else "MISMATCH"
        match = "MATCH" if our_db == vienna_db else "differ"
        print(f"{seq:<20} our={our_e:6.2f} vienna_eval_of_ours={vienna_e_of_our_structure:6.2f} "
              f"[{agree}]   our_db={our_db}  vienna_mfe_db={vienna_db} [{match}]")
