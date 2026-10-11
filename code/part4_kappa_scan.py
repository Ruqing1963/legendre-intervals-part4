"""
part4_kappa_scan.py -- Part IV, Section 3.2: the coset parameter kappa_0 of the classes 4 1^2 and 4 2 from ground-field data.

Parametrisation (Lemma 3.5): tau_0 = y^2 + y (y in F_q \\ {0,1}), rho_0 = s tau_0^2 (s in F_q^*), z = y^4 (so z^2 + z = tau_0^4),
b = rho_0 z / tau_0, a = b^2/rho_0 + tau_0^4 + rho_0 tau_0^2 + rho_0^2; then u^6 + a u^2 + b u + b^2 = (u^2 + tau_0 u + rho_0) P(u)
with the quartic cofactor P(u) = u^4 + tau_0 u^3 + (tau_0^2 + rho_0) u^2 + tau_0^3 u + b^2/rho_0.  Every pair (a, b) whose sextic
has a quadratic factor arises exactly once; the pair is of class 4 1^2 or 4 2 iff P is irreducible over F_q, and then the class
is 4 1^2 iff Tr(s) = Tr(rho_0/tau_0^2) = 0.
For each such (y, s) the script computes
  * kappa_0 by the trace expression (eq:q4r1) of Lemma 3.5 and by the closed form of Proposition 3.6,
        kappa_0 = Tr_{q/2}( z^3 + rho_0^2/z + (z+1) tau_0^2/rho_0 ),
  * the factorisation type of P and of the depressed quartic Phi(x) = x^4 + a' x^2 + b' x + c',  a' = 1/(1+t),
    b' = c' = 1/(s(1+t)), t = (y/(y+1))^4, obtained from P by u = tau_0 (x^{-1} + 1) (they must agree), and the resolvent
    criterion: P irreducible iff the cubic Y^3 + a' Y + b' has exactly one root r in F_q with Tr(c'/(a'+r^2)^2) = 1 and
    Tr(a'/r^2) = 0 (k odd),
and tabulates kappa_0 by class (and separately for the pure pairs a = 0, i.e. lambda_3 = 0).
Usage: python part4_kappa_scan.py 3 5 7
"""
import os
import sys
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from quintic_check import Poly  # noqa: E402
from legendre_d6_Wp_pattern import factor_degrees  # noqa: E402


def run(k):
    q = 1 << k
    P = Poly(q)
    A, M, INV = P.A, P.M, P.INV
    add = lambda a, b: int(A[a, b])
    mul = lambda a, b: int(M[a, b])
    inv = lambda a: int(INV[a])
    div = lambda a, b: mul(a, inv(b))

    def pw(a, e):
        r = 1
        for _ in range(e):
            r = mul(r, a)
        return r

    def sq(a):
        return mul(a, a)

    def tr(a):
        s, x = 0, a
        for _ in range(k):
            s = add(s, x)
            x = sq(x)
        assert s in (0, 1)
        return s

    def sqrt(a):
        return pw(a, q // 2) if a else 0

    by_class = Counter()
    checks = Counter()
    for y in range(2, q):
        tau0 = add(sq(y), y)
        e2abs = tr(add(y, pw(y, 3)))                       # e_2^abs(tau_0) = Tr(y + y^3)
        v = div(y, add(y, 1))
        t = pw(v, 4)
        for s in range(1, q):
            rho0 = mul(s, sq(tau0))
            z = pw(y, 4)
            b = div(mul(rho0, z), tau0)
            a = add(add(add(div(sq(b), rho0), pw(tau0, 4)), mul(rho0, sq(tau0))), sq(rho0))
            P4 = [div(sq(b), rho0), pw(tau0, 3), add(sq(tau0), rho0), tau0, 1]
            assert P.trim(P.mul(P4, [rho0, tau0, 1])) == P.trim([sq(b), b, a, 0, 0, 0, 1])
            typ = factor_degrees(P, P4, 4)
            # depressed quartic and its resolvent cubic
            d = mul(s, add(1, t))
            ap, bp = inv(add(1, t)), inv(d)
            Phi = [bp, bp, ap, 0, 1]
            typPhi = factor_degrees(P, Phi, 4)
            roots = [r for r in range(q) if add(add(pw(r, 3), mul(ap, r)), bp) == 0]
            crit = len(roots) == 1 and tr(div(bp, sq(add(ap, sq(roots[0]))))) == 1 and tr(div(ap, sq(roots[0]))) == 0
            checks[("Phi has the type of P", typ == typPhi)] += 1
            checks[("resolvent criterion <=> P irreducible", crit == (typ == (4,)))] += 1
            if typ != (4,):
                continue
            cls = "4 1^2" if tr(s) == 0 else "4 2"
            l5 = inv(b)
            l3 = mul(sqrt(a), l5)
            p3 = add(pw(tau0, 3), mul(rho0, tau0))
            p5 = add(add(pw(tau0, 5), mul(sq(rho0), tau0)), div(mul(sq(b), tau0), rho0))
            kappa_old = e2abs ^ tr(add(add(add(sq(tau0), rho0), mul(l5, p5)), mul(l3, p3)))
            kappa_new = tr(add(add(pw(z, 3), div(sq(rho0), z)), div(mul(add(z, 1), sq(tau0)), rho0)))
            checks[("closed form = eq:q4r1", kappa_old == kappa_new)] += 1
            by_class[(cls, "pure" if a == 0 else "mixed", "kappa0", kappa_old)] += 1
    print(f"\nk={k}, q={q}: (y,s) pairs with irreducible cofactor = pairs of the classes 4 1^2, 4 2")
    for key, n in sorted(checks.items(), key=str):
        print("   check", key, "x", n)
    tot = Counter()
    for (cls, pure, _, kap), n in sorted(by_class.items()):
        print(f"   class {cls:6s} {pure:5s} kappa0={kap}: {n}")
        tot[(cls, kap)] += n
    print("   totals by (class, kappa0):", dict(sorted(tot.items())))


if __name__ == "__main__":
    for k in [int(a) for a in sys.argv[1:]] or [3, 5]:
        run(k)
