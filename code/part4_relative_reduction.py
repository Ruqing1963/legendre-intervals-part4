"""
part4_relative_reduction.py -- Part IV, Lemma 3.x (reduction of s_4 to the quotient F_{q^4}/F_{q^2}).

On F_{q^2} the Z/4-form q^(4)_lambda is 2 Tr_{q^2/2} (multiplicity formula), i.e. twice a linear functional, and its
polar form is B(x,y) = 2 Tr_{q^4/2}(y T'_lambda(x)).  Summing over the cosets of F_{q^2} gives
    s_4(lambda) = sum_{c in A_lambda} i^{q^(4)_lambda(c)},    A_lambda = {c in F_{q^4} : Tr_{q^4/q^2}(T'_lambda(c)) = 1},
where A_lambda = T'^{-1}(F_{q^2} + t_1) (t_1 + t_1^{q^2} = 1) is a union of cosets of F_{q^2} on each of which q^(4) is
constant; A_lambda is a coset of T'^{-1}(F_{q^2}), of dimension 2k + 2 (K_lambda of dimension 4 plus 2k - 2 = dim of
T'(F_{q^4}) meet F_{q^2}).  Since |A_lambda| = 4q^2 and s_4 in {0, +-4q^2}: s_4 != 0 iff q^(4) is constant on A_lambda,
and then s_4 = 4q^2 i^{that constant} (the constant is 0 or 2).
The script checks the identity for every lambda_1 and every pair of the classes 4 1^2, 4 2 at q = 8 (and a sample at
q = 32), computing both sides over F_{2^{4k}}, and records dim A_lambda and the number of c in A_lambda with
q^(4)(c) = 0, 1, 2, 3.
Usage: python part4_relative_reduction.py 3          python part4_relative_reduction.py 5 4   (4 pairs per class)
"""
import os
import sys
from collections import Counter

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import legendre_var_check as LV  # noqa: E402
from quintic_check import Poly  # noqa: E402
from legendre_d6_Wp_pattern import factor_degrees  # noqa: E402
from part4_brown_invariants import embedding  # noqa: E402
import witt_lpoly  # noqa: E402  (installs the primitive polynomials)


def q4_values(k, l1, l3, l5, gf, emb):
    """array over x = 1..gf.size-1 of q^(4)_lambda(x) in Z/4 (same routine as witt_lpoly.witt_sums)."""
    n = 4 * k
    order = gf.order
    x = np.arange(1, gf.size, dtype=np.int64)
    logx = gf.log[x]

    def mono(c, e):
        if c == 0:
            return np.zeros_like(x)
        return gf.exp[(gf.log[c] + e * logx) % order]
    fx = mono(emb[l5], 5) ^ mono(emb[l3], 3) ^ mono(emb[l1], 1)
    trf = LV._abs_trace_parity(gf, fx.astype(np.uint32)).astype(np.int64)
    s0 = np.zeros_like(x); s1 = np.zeros_like(x)
    c = x.copy()
    for _ in range(n):
        prod = np.zeros_like(x)
        mask = (s0 != 0) & (c != 0)
        prod[mask] = gf.exp[(gf.log[s0[mask]] + gf.log[c[mask]]) % order]
        s1 ^= prod
        s0 ^= c
        c = gf.exp[(2 * gf.log[c]) % order]
    return (s0 + 2 * s1 + 2 * trf) % 4


def run(k, max_pairs=None):
    q = 1 << k
    P = Poly(q)
    gf = LV.GF(4 * k)
    emb = embedding(k, gf)
    order = gf.order
    x = np.arange(1, gf.size, dtype=np.int64)
    logx = gf.log[x]
    inv2 = pow(2, -1, order)
    inv4 = pow(4, -1, order)

    def mono(c, e):                      # c * x^e, e may be a modular inverse exponent (roots)
        if c == 0:
            return np.zeros_like(x)
        return gf.exp[(gf.log[c] + e * logx) % order]

    def frob_q2(v):
        out = np.zeros_like(v)
        m = v != 0
        out[m] = gf.exp[(q * q * gf.log[v[m]]) % order]
        return out

    i4 = np.array([1, 1j, -1, -1j])
    stats = Counter()
    done = Counter()
    for l5 in range(1, q):
        for l3 in range(q):
            B6 = [1, l5, int(P.M[l3, l3]), 0, 0, 0, int(P.M[l5, l5])]
            cls = factor_degrees(P, B6, 6)
            if cls not in ((1, 1, 4), (2, 4)):
                continue
            if max_pairs and done[cls] >= max_pairs:
                continue
            done[cls] += 1
            # T'(c) = l5 c^4 + l3 c^2 + c + (l3 c)^{1/2} + (l5 c)^{1/4}
            # roots inside F_q: l3^{1/2} = l3^{q/2}, l5^{1/4} = l5^{q/4}
            l3h = pow_in_Fq(P, l3, q // 2)
            l5qq = pow_in_Fq(P, l5, q // 4)       # l5^{1/4} = l5^{q/4}
            Tp = mono(emb[l5], 4) ^ mono(emb[l3], 2) ^ x ^ mono(emb[l3h], inv2) ^ mono(emb[l5qq], inv4)
            trT = Tp ^ frob_q2(Tp)                # Tr_{q^4/q^2} T'(c)
            A = trT == 1                          # the element 1 of gf is the integer 1
            dimA = int(round(np.log2(A.sum()))) if A.sum() else -1
            assert A.sum() == (1 << dimA)
            for l1 in range(q):
                vals = q4_values(k, l1, l3, l5, gf, emb)
                s4 = i4[vals].sum() + 1
                red = i4[vals[A]].sum()
                hist = tuple(int((vals[A] == j).sum()) for j in range(4))
                const = {j for j in range(4) if hist[j] == A.sum()}
                stats[(cls, "dim A", dimA, "s4 = sum_A i^q", abs(s4 - red) < 1e-6, "s4/4q^2", int(round(s4.real / (4 * q * q))),
                       "q^(4) constant on A", bool(const), "value", tuple(const), "counts of q^(4)=0,1,2,3 on A", hist)] += 1
    print(f"\nk={k}, q={q}: reduction of s_4 to the affine space A = {{c : Tr_(q^4/q^2) T'(c) = 1}}")
    for key, n in sorted(stats.items(), key=str):
        print("  ", key, "x", n)


def pow_in_Fq(P, a, e):
    r = 1
    for _ in range(e):
        r = int(P.M[r, a])
    return r


if __name__ == "__main__":
    k = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    mp = int(sys.argv[2]) if len(sys.argv) > 2 else None
    run(k, mp)
