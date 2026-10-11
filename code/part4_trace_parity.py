"""
part4_trace_parity.py -- Part IV, Theorem 3.x (trace parity on the coset T'^{-1}(1)), the case m = 1.

Claim (k odd): every rational root r of B_lambda (equivalently of u^6 + a u^2 + b u + b^2) satisfies
    Tr_{q/2}(r) = f_1(C) + 1 (mod 2),   f_1(C) = dim K_lambda^sigma = log_2(1 + #duads fixed by pi_lambda),
i.e. Tr(r) = 0 in the classes 2 1^4 (f_1 = 3), 3 2 1 (f_1 = 1), 4 1^2 (f_1 = 1) and Tr(r) = 1 in the classes
1^6 (f_1 = 4), 2^2 1^2 (f_1 = 2), 3 1^3 (f_1 = 2), 5 1 (f_1 = 0).  (General statement: Tr_{q^m/2}(w) = mk + f_m mod 2 for
every w in T'^{-1}(1) meeting F_{q^m}, provided s_m(lambda_1) != 0 for some lambda_1.)
The script lists, for q = 8, 32, 128 and every (a, b), the traces of the rational roots by factorisation type.
Usage: python part4_trace_parity.py 3 5 7
"""
import os
import sys
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from quintic_check import Poly  # noqa: E402
from legendre_d6_Wp_pattern import factor_degrees  # noqa: E402

F1 = {(1, 1, 1, 1, 1, 1): 4, (1, 1, 1, 1, 2): 3, (1, 1, 2, 2): 2, (1, 1, 1, 3): 2, (1, 2, 3): 1, (1, 1, 4): 1, (1, 5): 0}


def run(k):
    q = 1 << k
    P = Poly(q)
    A, M, INV = P.A, P.M, P.INV

    def tr(a):
        s, x = 0, a
        for _ in range(k):
            s = int(A[s, x]); x = int(M[x, x])
        return s

    # rational roots: for each (u, b) the unique a with F_{a,b}(u) = 0 is a = u^4 + b/u + b^2/u^2
    roots = {}
    for b in range(1, q):
        for u in range(1, q):
            iu = int(INV[u])
            a = int(A[int(A[int(M[int(M[u, u]), int(M[u, u])]), int(M[b, iu])]), int(M[int(M[b, b]), int(M[iu, iu])])])
            roots.setdefault((a, b), []).append(u)
    stats = Counter()
    for (a, b), us in roots.items():
        F6 = [int(M[b, b]), b, a, 0, 0, 0, 1]
        cls = factor_degrees(P, F6, 6)
        assert len(us) == sum(1 for d in cls if d == 1)
        for u in us:
            stats[(cls, "f_1", F1[cls], "Tr(r)", tr(u))] += 1
    print(f"\nk={k}, q={q}: traces of the rational roots by class (claim: Tr(r) = k + f_1 mod 2, i.e. f_1 + {k % 2 ^ 1} here)")
    for key, n in sorted(stats.items(), key=str):
        print("  ", key, "x", n, "  ok" if key[4] == (key[2] + k) % 2 else "  MISMATCH")


if __name__ == "__main__":
    for k in [int(a) for a in sys.argv[1:]] or [3, 5]:
        run(k)
