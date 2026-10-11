"""
part4_sign_rules.py -- Part IV, Theorems 3.x (signs of s_4 for N_0 | 2 and of s_6 for N_0 in {3,6}) from the hexad formula.

For a pair (lambda_3, lambda_5) with N_0 | 2 (classes 1^6, 2 1^4, 2^2 1^2, 2^3) put phi_i = q^(2)_lambda(r_i)/2 in F_2 for the six
roots (all q^(2)(r_i) are even), n_1 = #{i : phi_i = 1}.  Then (Theorem): sum phi_i = 0, the value at a total r_i+r_j+r_l is
phi_i+phi_j+phi_l, and
    s_4 = q^2 ( sum_i (-1)^{phi_i} - sum_{totals} (-1)^{phi} ) = +4q^2 if n_1 = 2, 6  and  -4q^2 if n_1 = 0, 4.
Rational roots have phi = Tr_{q/2}(r) = f_1 + 1 (trace parity), conjugate roots have equal phi; this gives
  1^6: n_1 = 6 (+);  2^3: an odd number of pairs has phi = 1 because the rational totals have phi = 1 (+);
  2 1^4: sign = (-1)^{1 + phi(pair)};  2^2 1^2: sign = (-1)^{phi(pair 1) + phi(pair 2)}.
For N_0 in {3,6} (classes 3 1^3, 3^2, 6; coset points in F_{q^3} have odd Witt value) put phi(w) = (q^(3)(w) - 1)/2; then
s_6 = i q^3 ( sum_roots (-1)^phi - sum_totals (-1)^phi ) and the same combinatorics gives s_6 = +-4 i q^3, hence S'' = +4q^6,
provided s_3(lambda_1) != 0 for some lambda_1 (hypothesis H; automatic for 3 1^3 and 3^2, checked here for 6).
The script computes, for every pair at q = 8 (and a sample at q = 32): the phi's with the Frobenius orbit structure,
n_1, the predicted and the actual s_{2m} (vectorised Witt sums), and H (= q^(3)_{lambda_1=0} vanishes on K cap F_{q^3}).
Usage: python part4_sign_rules.py 3          python part4_sign_rules.py 5 3    (q=32, 3 pairs per class, m=2 only)
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

LV.PRIMITIVE.setdefault(18, (1 << 18) | (1 << 7) | 1)
N0 = {(1, 1, 1, 1, 1, 1): 1, (1, 1, 1, 1, 2): 2, (1, 1, 2, 2): 2, (2, 2, 2): 2, (1, 1, 1, 3): 3, (3, 3): 3, (1, 2, 3): 6,
      (1, 1, 4): 4, (2, 4): 4, (1, 5): 5, (6,): 6}


def run(k, max_pairs=None):
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
    ms = (2, 3) if k == 3 else (2,)
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

    def orbit_len(x):
        y, n = Fd.frob(x), 1
        while not all(a == b for a, b in zip(y, x)):
            y = Fd.frob(y); n += 1
        return n

    stats = Counter()
    done = Counter()
    for l5 in range(1, q):
        for l3 in range(q):
            B6 = [1, l5, int(P.M[l3, l3]), 0, 0, 0, int(P.M[l5, l5])]
            cls = factor_degrees(P, B6, 6)
            n0 = N0[cls]
            if n0 in (4, 5) or cls == (1, 2, 3):
                continue
            m = 2 if n0 in (1, 2) else 3
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
            roots = [w for w in coset if all(v == 0 for v in Bval(w))]
            totals = [w for w in coset if not all(v == 0 for v in Bval(w))]
            assert len(roots) == 6 and len(totals) == 10
            Kel = [Fd.K12.frombits(sum((int(c) * b for c, b in zip(cf, K)), np.zeros(N, dtype=np.uint8)) % 2)
                   for cf in itertools.product((0, 1), repeat=4)]
            # hypothesis H for m = 3: q^(3)_0 vanishes on K cap F_{q^3}
            H = None
            if m == 3:
                H = all(qm(kappa, 3, 0, l3, l5) == 0 for kappa in Kel if orbit_len(kappa) in (1, 3))
            roots_m = [w for w in roots if m % orbit_len(w) == 0]
            totals_m = [w for w in totals if m % orbit_len(w) == 0]
            for l1 in range(q):
                vals_r = [qm(w, m, l1, l3, l5) for w in roots_m]
                vals_t = [qm(w, m, l1, l3, l5) for w in totals_m]
                parity_ok = all(v % 2 == (0 if m == 2 else 1) for v in vals_r + vals_t)
                phi_r = [(v - (v % 2)) // 2 % 2 for v in vals_r]
                phi_t = [(v - (v % 2)) // 2 % 2 for v in vals_t]
                n1 = sum(phi_r)
                orbit_phi = tuple(sorted((orbit_len(w), ph) for w, ph in zip(roots_m, phi_r)))
                if len(roots_m) == 6:                                # sixteen points: the n_1 rule
                    pred = 4 if n1 % 4 == 2 else -4
                else:                                                # class 6: four totals, no roots
                    pred = -4 * (-1) ** phi_t[0] if len(set(phi_t)) == 1 else 0
                unit = 1 if m == 2 else 1j
                lhs = witt_sums(k, 2 * m, l1, l3, l5, gfs[m], embs[m])
                actual = lhs / q ** m
                ok = abs(actual - pred * unit) < 1e-6
                key = (cls, "m", m, "parity ok", parity_ok, "sum phi = 0", n1 % 2 == 0, "n1", n1,
                       "pred s_2m/q^m", pred * unit, "actual", complex(round(actual.real), round(actual.imag)), "rule ok", ok)
                if m == 3:
                    key = key + ("H: q3 = 0 on K cap F_q3", H)
                if q <= 8:
                    key = key + ("(orbit length, phi) of roots", orbit_phi)
                stats[key] += 1
    print(f"\nk={k}, q={q}: sign rules from the hexad formula (N_0 | 2 with m = 2; N_0 in {{3,6}} except 3 2 1 with m = 3)")
    for key, n in sorted(stats.items(), key=str):
        print("  ", key, "x", n)


if __name__ == "__main__":
    k = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    mp = int(sys.argv[2]) if len(sys.argv) > 2 else None
    run(k, mp)
