"""
part4_hexad_formula.py -- Part IV, Theorem 3.x (the hexad formula for the small-field Gauss sums s_{2m}).

For lambda_5 != 0 let T'^{-1}(1) = r + K_lambda be the coset of the kernel containing the six roots of B_lambda (the other
ten elements are the synthematic totals).  Summing the Z/4-form q^(2m)_lambda over the cosets of F_{q^m} in F_{q^{2m}}
(on F_{q^m} the form is 2 Tr_{q^m/2}) and evaluating it at the solutions c of c + c^{q^m} = w gives
    s_{2m}(lambda) = q^m * sum_{w in T'^{-1}(1), w in F_{q^m}}  eps(w) * i^{q^(m)_lambda(w)},
    eps(w) = +1 if B_lambda(w) = 0 (w a root), -1 if w is a synthematic total,
with q^(m)_lambda(w) = T_{mk}(w) + 2 Tr_{q^m/2}(lambda_5 w^5 + lambda_3 w^3 + lambda_1 w).
The script checks this for m = 2 (all pairs with N_0 | 4 at q = 8; a sample per class at q = 32) and for m = 3 (all pairs
with N_0 in {3,6} at q = 8), for every lambda_1: the left side by vectorised Witt sums over F_{2^{2mk}}, the right side from
the hexad in F_{q^12}.  It also records the multiset of the signed terms.
Usage: python part4_hexad_formula.py 3            python part4_hexad_formula.py 5 2 3   (k=5, m=2 only, 3 pairs per class)
"""
import itertools
import os
import sys
from collections import Counter

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import legendre_d6_qf as QF  # noqa: E402
import legendre_var_check as LV  # noqa: E402
from hexad_search import Field  # noqa: E402
from quintic_check import Poly  # noqa: E402
from legendre_d6_Wp_pattern import factor_degrees  # noqa: E402
from part4_brown_invariants import embedding  # noqa: E402
from witt_lpoly import witt_sums  # noqa: E402

LV.PRIMITIVE.setdefault(18, (1 << 18) | (1 << 7) | 1)          # x^18 + x^7 + 1 (primitive; checked by the constructor)

N0 = {(1, 1, 1, 1, 1, 1): 1, (1, 1, 1, 1, 2): 2, (1, 1, 2, 2): 2, (2, 2, 2): 2, (1, 1, 1, 3): 3, (3, 3): 3, (1, 2, 3): 6,
      (1, 1, 4): 4, (2, 4): 4, (1, 5): 5, (6,): 6}


def run(k, ms=(2, 3), max_pairs=None):
    Fd = Field(k)
    q = Fd.q
    F = Fd.F
    P = Poly(q)
    SQ, SQRT = Fd.SQ, Fd.SQRT
    SQ2 = SQ.dot(SQ) % 2
    SQRT2 = SQRT.dot(SQRT) % 2
    N = Fd.N
    I = np.eye(N, dtype=np.uint8)
    cache = {}
    gfs = {m: LV.GF(2 * m * k) for m in ms}
    embs = {m: embedding(k, gfs[m]) for m in ms}

    def S(c):
        if c not in cache:
            cache[c] = Fd.K12.linmap_matrix(lambda a: Fd.scal(c, a))
        return cache[c]

    def in_Fq_val(a):
        assert Fd.in_Fq(a), a
        return int(a[0])

    def witt_T(x, nconj):
        conj = [x]
        for _ in range(nconj - 1):
            conj.append(Fd.mul(conj[-1], conj[-1]))
        s0, s1 = Fd.zero, Fd.zero
        for c in conj:
            s1 = Fd.add(s1, Fd.mul(s0, c))
            s0 = Fd.add(s0, c)
        return (in_Fq_val(s0) + 2 * in_Fq_val(s1)) % 4

    def tr_rel(x, m):
        t, y = Fd.zero, x
        for _ in range(m):
            t = Fd.add(t, y); y = Fd.frob(y)
        return in_Fq_val(t)

    def qm(x, m, l1, l3, l5):
        x2 = Fd.mul(x, x); x3 = Fd.mul(x2, x); x5 = Fd.mul(x3, x2)
        f = Fd.add(Fd.add(Fd.scal(l5, x5), Fd.scal(l3, x3)), Fd.scal(l1, x))
        return (witt_T(x, m * k) + 2 * F.tr2(tr_rel(f, m))) % 4

    def in_subfield(x, m):
        y = x
        for _ in range(m):
            y = Fd.frob(y)
        return y == x if isinstance(y, int) else all(a == b for a, b in zip(y, x))

    stats = Counter()
    done = Counter()
    for l5 in range(1, q):
        for l3 in range(q):
            B6 = [1, l5, int(P.M[l3, l3]), 0, 0, 0, int(P.M[l5, l5])]
            cls = factor_degrees(P, B6, 6)
            n0 = N0[cls]
            if n0 == 5:                       # class 5 1: the hexad lies in F_{q^5}, outside F_{q^12}; S'' is proved there
                continue
            m = 2 if n0 in (1, 2, 4) else 3
            if m not in ms:
                continue
            if max_pairs and done[cls] >= max_pairs:
                continue
            done[cls] += 1
            Tp = (S(l5).dot(SQ2) + S(l3).dot(SQ) + I + S(int(F.SQRT[l3])).dot(SQRT)
                  + S(int(F.SQRT[int(F.SQRT[l5])])).dot(SQRT2)) % 2
            K = QF.f2_nullspace(SQ2.dot(Tp) % 2)
            one = Fd.K12.tobits(Fd.one)
            sol = QF.f2_nullspace(np.concatenate([Tp, one[:, None]], axis=1))
            x0 = next(v[:-1] for v in sol if v[-1] == 1)
            coset = []
            for cf in itertools.product((0, 1), repeat=4):
                v = x0.copy()
                for c, b in zip(cf, K):
                    if c:
                        v ^= b
                coset.append(Fd.K12.frombits(v % 2))
            sq3, sq5 = int(Fd.M[l3, l3]), int(Fd.M[l5, l5])

            def Bval(u):
                u2 = Fd.mul(u, u); u6 = Fd.mul(Fd.mul(u2, u2), u2)
                return Fd.add(Fd.add(Fd.scal(sq5, u6), Fd.scal(sq3, u2)), Fd.add(Fd.scal(l5, u), Fd.one))
            W = [(w, all(v == 0 for v in Bval(w))) for w in coset if in_subfield(w, m)]
            nroots = sum(1 for _, isroot in W if isroot)
            ok_all = True
            terms_all = Counter()
            for l1 in range(q):
                lhs = witt_sums(k, 2 * m, l1, l3, l5, gfs[m], embs[m])
                terms = [((1 if isroot else -1), qm(w, m, l1, l3, l5)) for w, isroot in W]
                rhs = q ** m * sum(e * (1j ** v) for e, v in terms)
                ok = abs(lhs - rhs) < 1e-6
                ok_all &= ok
                terms_all[(int(round(lhs.real / q ** m)), int(round(lhs.imag / q ** m)), tuple(sorted(terms)))] += 1
            stats[(cls, "m", m, "#T'^{-1}(1) in F_{q^m}", len(W), "roots", nroots, "formula ok for all l1", ok_all)] += 1
            if q <= 8:
                for key, n in terms_all.items():
                    stats[(cls, "m", m, "s_{2m}/q^m", key[0] + 1j * key[1], "signed terms (eps, q^(m)(w))", key[2])] += n
    print(f"\nk={k}, q={q}: hexad formula s_{{2m}} = q^m sum_{{w in T'^(-1)(1) cap F_(q^m)}} eps(w) i^{{q^(m)(w)}}")
    for key, n in sorted(stats.items(), key=str):
        print("  ", key, "x", n)


if __name__ == "__main__":
    k = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    ms = tuple(int(a) for a in sys.argv[2:-1]) if len(sys.argv) > 3 else ((2, 3) if len(sys.argv) <= 2 else (int(sys.argv[2]),))
    mp = int(sys.argv[-1]) if len(sys.argv) > 3 else None
    run(k, ms, mp)
