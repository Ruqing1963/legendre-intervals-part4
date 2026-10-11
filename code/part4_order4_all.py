"""
part4_order4_all.py -- Part IV, Section 3.2: every pair of the classes 4 1^2 and 4 2 at q = 2^k.

For each (lambda_3, lambda_5) with B_lambda of class 4 1^2 or 4 2 the script finds the hexad in F_{q^12}, the 4-cycle
root r_1 and the pair {r_5, r_6} (tau_0 = r_5 + r_6, rho_0 = r_5 r_6, b = 1/lambda_5), and records
  * I3 = q^(4)_{lambda_1=0}(r_1) in Z/4 (direct Witt computation) and kappa_0 = q^(4)(r_1 + r_5)/2 = I3/2;
  * the trace expression (eq:q4r1) of Lemma 3.5 for I3/2, and the reduced closed form of Proposition 3.6,
        kappa_0 = Tr_{q/2}( z^3 + rho_0^2/z + (z+1) tau_0^2/rho_0 ),   z = b tau_0/rho_0  (z^2 + z = tau_0^4);
  * (unless --fast) s_4(lambda_1) for all lambda_1 by vectorised Witt sums over F_{2^{4k}}: the set where s_4 != 0
    (checked to be the affine hyperplane Tr(lambda_1 tau_0) = kappa_0), the sign of s_4 there, the comparison with
    (-1)^{Tr(rho_0/tau_0^2)}, and W' = q^4 sum_{lambda_1} s_4.
Both values kappa_0 = 0 and kappa_0 = 1 occur (q = 32: 110 and 130 pairs); W' = +-2q^7 in either case.
Usage: python part4_order4_all.py 3            python part4_order4_all.py 5           python part4_order4_all.py 5 --fast
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


def run(k, fast=False):
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
    gf4 = emb = None
    if not fast:
        gf4 = LV.GF(4 * k)
        emb = embedding(k, gf4)

    def S(c):
        if c not in cache:
            cache[c] = Fd.K12.linmap_matrix(lambda a: Fd.scal(c, a))
        return cache[c]

    def in_Fq_val(a):
        assert Fd.in_Fq(a), a
        return int(a[0])

    def witt_Z4(x, l1, l3, l5):
        """q^(4)(x) for x in F_{q^4} (given inside F_{q^12}) via q^(12)(x) = q^(4)(x) + 2 Tr_{q^4/2}(x)."""
        conj = [x]
        for _ in range(12 * k - 1):
            conj.append(Fd.mul(conj[-1], conj[-1]))
        s0, s1 = Fd.zero, Fd.zero
        for c in conj:
            s1 = Fd.add(s1, Fd.mul(s0, c))
            s0 = Fd.add(s0, c)
        T = (in_Fq_val(s0) + 2 * in_Fq_val(s1)) % 4
        x2 = Fd.mul(x, x); x3 = Fd.mul(x2, x); x5 = Fd.mul(x3, x2)
        f = Fd.add(Fd.add(Fd.scal(l5, x5), Fd.scal(l3, x3)), Fd.scal(l1, x))
        q12 = (T + 2 * Fd.abs_trace(f)) % 4
        t, y = Fd.zero, x
        for _ in range(4):
            t = Fd.add(t, y)
            y = Fd.frob(y)
        return (q12 - 2 * F.tr2(in_Fq_val(t))) % 4

    fm = lambda a, b: int(F.M[a, b])
    finv = lambda a: int(F.INV[a])

    def e2abs(t0):
        conj, x = [], t0
        for _ in range(k):
            conj.append(x); x = int(F.SQ[x])
        e = 0
        for i in range(k):
            for j in range(i + 1, k):
                e ^= fm(conj[i], conj[j])
        assert e in (0, 1)
        return e

    stats = Counter()
    kappa_by_class = Counter()
    rows = []
    for l5 in range(1, q):
        for l3 in range(q):
            B6 = [1, l5, fm(l3, l3), 0, 0, 0, fm(l5, l5)]
            cls = factor_degrees(P, B6, 6)
            if cls not in ((1, 1, 4), (2, 4)):
                continue
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
            sq3, sq5 = fm(l3, l3), fm(l5, l5)

            def Bval(u):
                u2 = Fd.mul(u, u); u6 = Fd.mul(Fd.mul(u2, u2), u2)
                return Fd.add(Fd.add(Fd.scal(sq5, u6), Fd.scal(sq3, u2)), Fd.add(Fd.scal(l5, u), Fd.one))
            roots = [u for u in coset if all(v == 0 for v in Bval(u))]
            assert len(roots) == 6
            orb, seen = {}, set()
            for r in roots:
                t = tuple(r)
                if t in seen:
                    continue
                o, y = [], t
                while y not in seen:
                    seen.add(y); o.append(y); y = tuple(Fd.frob(list(y)))
                for z in o:
                    orb[z] = len(o)
            four = [list(r) for r in roots if orb[tuple(r)] == 4]
            other = [list(r) for r in roots if orb[tuple(r)] != 4]
            r1, r5, r6 = four[0], other[0], other[1]
            t0 = in_Fq_val(Fd.add(r5, r6))
            p = in_Fq_val(Fd.mul(r5, r6))
            I1 = F.tr2(fm(p, finv(fm(t0, t0))))
            I3 = witt_Z4(r1, 0, l3, l5)
            kappa0 = witt_Z4(Fd.add(r1, r5), 0, l3, l5)
            assert kappa0 % 2 == 0 and kappa0 == I3
            kappa0 //= 2
            b_ = finv(l5)
            t2 = fm(t0, t0); t3 = fm(t2, t0); t5 = fm(t3, t2)
            p3 = t3 ^ fm(p, t0)
            b2 = fm(b_, b_)
            p5 = t5 ^ fm(fm(p, p), t0) ^ fm(fm(b2, t0), finv(p))
            e2 = e2abs(t0)
            formula_old = (e2 + F.tr2(t2 ^ p ^ fm(l5, p5) ^ fm(l3, p3))) % 2
            z = fm(fm(b_, t0), finv(p))                       # z = b tau0 / rho0,  z^2 + z = tau0^4
            assert fm(z, z) ^ z == fm(t2, t2)
            formula_new = F.tr2(fm(fm(z, z), z) ^ fm(fm(p, p), finv(z)) ^ fm(fm(z ^ 1, t2), finv(p)))
            rec = (cls, "kappa0", kappa0, "eq:q4r1 ok", formula_old == kappa0, "closed form ok", formula_new == kappa0,
                   "I1=Tr(rho0/tau0^2)", I1)
            kappa_by_class[(cls, kappa0)] += 1
            if not fast:
                s4 = {l1: witt_sums(k, 4, l1, l3, l5, gf4, emb) for l1 in range(q)}
                nz = {l1: v for l1, v in s4.items() if abs(v) > 1e-9}
                assert all(abs(v.imag) < 1e-9 and abs(abs(v.real) - 4 * q * q) < 1e-6 for v in nz.values())
                signs = {int(round(v.real / (4 * q * q))) for v in nz.values()}
                coset_ok = all((F.tr2(fm(l1, t0)) == kappa0) == (l1 in nz) for l1 in range(q))
                Wp = sum(s4.values()) * q ** 4
                rec = rec + ("#nonzero", len(nz), "nonzero iff Tr(l1 tau0)=kappa0", coset_ok, "sign s_4", tuple(sorted(signs)),
                             "sign=(-1)^I1", signs == {(-1) ** I1}, "W'/q^7", int(round(Wp.real / q ** 7)))
                rows.append((l3, l5, cls, kappa0, I1, tuple(sorted(signs)), int(round(Wp.real / q ** 7))))
            stats[rec] += 1
    print(f"\nk={k}, q={q}: classes 4 1^2 = (1,1,4) and 4 2 = (2,4){' (fast: no s_4)' if fast else ''}")
    for key, n in sorted(stats.items(), key=str):
        print(f"   {key}  x{n}")
    print("   kappa_0 by class:", dict(sorted(kappa_by_class.items())))
    if rows:
        print("   per pair (l3, l5, class, kappa0, I1, signs, W'/q^7):")
        for r in rows:
            print("     ", r)


if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    run(int(args[0]) if args else 3, fast="--fast" in sys.argv)
