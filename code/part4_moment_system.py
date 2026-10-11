"""
part4_moment_system.py -- Part IV, Section 2.4: the elimination system of the moments.

With n_C(q) the class counts of Part III, Table 2 (mixed pairs lambda_3 != 0 and pure pairs lambda_3 = 0) and
f_C^(m) the number of fixed points of pi^m for a permutation pi of class C in S_6, the m-th moment of the sextic family is
    sum_{(a,b)} N_m(a,b) = sum_C n_C(q) f_C^(m).
The script prints these polynomials for m = 1..6, the implied values of the Gauss-sum totals
    Sigma_m(q) = sum_{theta in Theta_m} G_{H_m}(theta) = q^{m-2} ( sum N_m - (q-1)^2 + q [2|m] )
of Proposition 2.4, the rank of the linear system {f_C^(m)}_{m=1..6} in the eleven unknowns n_C (the pure counts are
theorems), and the quadratic interpolation through q = 8, 32, 128 of the per-type counts of Proposition 2.2(iv) and of the
kappa_0 counts of Proposition 3.6 (to show which of them are, and which are not, polynomial with simple coefficients).
Usage: python part4_moment_system.py
"""
from fractions import Fraction

import sympy as sp

q = sp.symbols("q")
CLASSES = ["1^6", "2 1^4", "2^2 1^2", "2^3", "3 1^3", "3^2", "3 2 1", "4 1^2", "4 2", "5 1", "6"]
MIXED = [(q - 2) * (q - 8) / 720, (q - 2) * (q - 8) / 48, q * (q - 2) / 16, (q - 2) * (q - 8) / 48,
         (q + 1) * (q - 2) / 18, (q + 1) * (q - 2) / 18, (q + 1) * (q - 2) / 6, (q - 2) * (q - 4) / 8, q * (q - 2) / 8,
         (q ** 2 + 1) / 5, (q + 1) * (q - 2) / 6]
PURE = [0, 0, 0, (q - 2) / 6, 0, 0, 0, (q - 2) / 2, 0, 0, (q + 1) / 3]
CYCLES = [(1, 1, 1, 1, 1, 1), (2, 1, 1, 1, 1), (2, 2, 1, 1), (2, 2, 2), (3, 1, 1, 1), (3, 3), (3, 2, 1), (4, 1, 1), (4, 2),
          (5, 1), (6,)]


def fixed(cycle_type, m):
    return sum(c for c in cycle_type if m % c == 0)


def main():
    F = sp.Matrix([[fixed(ct, m) for ct in CYCLES] for m in range(1, 7)])
    print("fixed-point matrix f_C^(m), rows m = 1..6, columns", CLASSES)
    sp.pprint(F)
    print("rank over Q:", F.rank(), " (eleven unknowns: the six moments alone leave", 11 - F.rank(), "degrees of freedom)")
    print("\npredicted moments sum N_m and Gauss-sum totals Sigma_m:")
    for m in range(1, 7):
        tot = sp.expand(sum((MIXED[i] + PURE[i]) * fixed(CYCLES[i], m) for i in range(11)))
        sigma = sp.expand(q ** (m - 2) * (tot - (q - 1) ** 2 + q * (1 if m % 2 == 0 else 0)))
        print(f"  m={m}:  sum N_m = {sp.factor(tot)} = {tot};   Sigma_m = {sp.factor(sigma)}")
        for qq in (8, 32, 128):
            pass
    print("\nnumerical values at q = 8, 32, 128:")
    for m in range(1, 7):
        tot = sum((MIXED[i] + PURE[i]) * fixed(CYCLES[i], m) for i in range(11))
        print(f"  m={m}: ", [int(tot.subs(q, qq)) for qq in (8, 32, 128)])

    def interp(vals, name):
        pts = [(8, vals[0]), (32, vals[1]), (128, vals[2])]
        poly = sp.interpolate(pts, q)
        print(f"  {name}: values {vals} -> quadratic through them: {sp.factor(poly)}  =  {sp.expand(poly)}")

    print("\nquadratic interpolation of the per-Gold-type counts of theta in Theta (moment3_structure.py):")
    interp([54, 540, 8442], "type (1,4)   ")
    interp([9, 333, 5337], "type (2,3)   ")
    interp([0, 150, 2604], "type (1,1,1,2)")
    interp([27, 495, 8127], "G_H != 0 (all types)")
    print("\nquadratic interpolation of the kappa_0 counts (part4_kappa_scan.py):")
    interp([6, 40, 924], "4 1^2, kappa_0 = 0")
    interp([0, 80, 1092], "4 1^2, kappa_0 = 1")
    interp([6, 70, 1106], "4 2,   kappa_0 = 0")
    interp([0, 50, 910], "4 2,   kappa_0 = 1")
    interp([3, 5, 21], "pure 4 1^2, kappa_0 = 0")
    interp([0, 10, 42], "pure 4 1^2, kappa_0 = 1")


if __name__ == "__main__":
    main()
