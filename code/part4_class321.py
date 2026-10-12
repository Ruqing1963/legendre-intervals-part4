"""
part4_class321.py -- Part IV, Theorem 3.21: the sign law in the class 3 2 1.

Claim.  Let pi_lambda have type 3 2 1, r4, r5 the two conjugate roots of B_lambda forming the 2-cycle, t0 = r4 + r5.
Then t0 is the unique nonzero element of ker T'_lambda cap F_q, and
      S''(lambda_1) = A_lambda * (-1)^{Tr_{q/2}(lambda_1 t0)},      A_lambda = +-4q^6 independent of lambda_1,
so S'' = +4q^6 for exactly q/2 values of lambda_1 and W' = 0.

Independent check: S'' = s_12 is computed from the L-polynomial of the rank-4 Witt sheaf, i.e. from the Witt sums
s_1..s_4 over F_{q^m}, m <= 4 (witt_lpoly.witt_sums, newton_p12) -- not from the hexad formula used in the proof.
t0 is found by brute force as the nonzero x in F_q with T'(x) = 0.
The class 6 is checked in the same way: S''(lambda_1) = +4q^6 for every lambda_1 (no hypothesis (H)).
Usage: python part4_class321.py 3 5      (q = 8: all pairs; q = 32: the first 12 pairs of each class)
"""
import os
import sys
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import legendre_d6_qf as QF  # noqa: E402
import legendre_var_check as LV  # noqa: E402
from quintic_check import Poly  # noqa: E402
from legendre_d6_Wp_pattern import factor_degrees  # noqa: E402
from witt_lpoly import witt_sums, newton_p12  # noqa: E402


def embeddings(k, gfs):
    """F_q (fqlib encoding) -> each gf, as in witt_lpoly.run"""
    q = 1 << k
    F = QF.Fq(q)
    emb = {}
    for m, gf in gfs.items():
        step = gf.order // (q - 1)
        target = None
        for t in range(q - 1):
            cand = int(gf.exp[(t * step) % gf.order])
            tab, a_pow, c_pow, ok = {}, 1, 1, True
            for _ in range(q - 1):
                tab[a_pow] = c_pow
                a_pow = int(F.M[a_pow, 2])
                c_pow = int(gf.exp[(gf.log[c_pow] + gf.log[cand]) % gf.order]) if c_pow else 0
            tab[0] = 0
            for a in range(1, min(q, 64)):
                for b in range(1, min(q, 64)):
                    if tab[int(F.A[a, b])] != (tab[a] ^ tab[b]):
                        ok = False
                        break
                if not ok:
                    break
            if ok:
                target = tab
                break
        emb[m] = target
    return emb


def run_class6(k, max_pairs=None):
    """class 6: S''(lambda_1) = +4q^6 for every lambda_1 (Theorem 3.21(ii)); independent L-polynomial check"""
    q = 1 << k
    F = QF.Fq(q)
    P = Poly(q)
    gfs = {m: LV.GF(m * k) for m in (1, 2, 3, 4)}
    emb = embeddings(k, gfs)
    unit = 4 * q ** 6
    stats = Counter()
    npairs = 0
    for l5 in range(1, q):
        for l3 in range(q):
            B6 = [1, l5, int(P.M[l3, l3]), 0, 0, 0, int(P.M[l5, l5])]
            if factor_degrees(P, B6, 6) != (6,):
                continue
            if max_pairs and npairs >= max_pairs:
                continue
            npairs += 1
            vals = set()
            for l1 in range(q):
                s = {m: witt_sums(k, m, l1, l3, l5, gfs[m], emb[m]) for m in (1, 2, 3, 4)}
                p12, _ = newton_p12(s)
                v = -p12 / unit
                vals.add(complex(round(v.real, 6), round(v.imag, 6)))
            stats[("S''/4q^6 over all l1", tuple(sorted(vals, key=str)))] += 1
    print(f"k={k}, q={q}: class 6 (mixed pairs), {npairs} pairs, every lambda_1")
    for key, n in sorted(stats.items(), key=str):
        print("  ", key, "x", n)


def run(k, max_pairs=None):
    q = 1 << k
    F = QF.Fq(q)
    P = Poly(q)
    gfs = {m: LV.GF(m * k) for m in (1, 2, 3, 4)}
    emb = embeddings(k, gfs)
    unit = 4 * q ** 6

    def Tp(x, l3, l5):
        x2 = int(F.M[x, x]); x4 = int(F.M[x2, x2])
        v = int(F.M[l5, x4]) ^ int(F.M[l3, x2]) ^ x
        v ^= int(F.SQRT[int(F.M[l3, x])])
        v ^= int(F.SQRT[int(F.SQRT[int(F.M[l5, x])])])
        return v

    stats = Counter()
    npairs = 0
    for l5 in range(1, q):
        for l3 in range(q):
            B6 = [1, l5, int(P.M[l3, l3]), 0, 0, 0, int(P.M[l5, l5])]
            if factor_degrees(P, B6, 6) != (1, 2, 3):
                continue
            if max_pairs and npairs >= max_pairs:
                continue
            npairs += 1
            ker = [x for x in range(1, q) if Tp(x, l3, l5) == 0]
            assert len(ker) == 1, ker
            t0 = ker[0]
            vals = []
            for l1 in range(q):
                s = {m: witt_sums(k, m, l1, l3, l5, gfs[m], emb[m]) for m in (1, 2, 3, 4)}
                p12, _ = newton_p12(s)
                vals.append(-p12)
            norm = [complex(round(v.real / unit, 6), round(v.imag / unit, 6)) for v in vals]
            chi = [(-1) ** F.tr2(int(F.M[l1, t0])) for l1 in range(q)]
            A = {norm[l1] * chi[l1] for l1 in range(q)}
            ok = len(A) == 1 and next(iter(A)) in (1, -1)
            nplus = sum(1 for v in norm if v == 1)
            stats[("law S''(l1) = A (-1)^Tr(l1 t0) holds", ok, "A", next(iter(A)) if len(A) == 1 else A,
                   "#{S''=+4q^6}", nplus, "of", q, "Tr(t0)", F.tr2(t0))] += 1
    print(f"k={k}, q={q}: class 3 2 1, {npairs} pairs, every lambda_1")
    for key, n in sorted(stats.items(), key=str):
        print("  ", key, "x", n)


if __name__ == "__main__":
    args = [int(a) for a in sys.argv[1:]] or [3, 5]
    for k in args:
        mp = None if k <= 3 else (12 if k == 5 else 3)
        run(k, max_pairs=mp)
        run_class6(k, max_pairs=mp)
        sys.stdout.flush()
