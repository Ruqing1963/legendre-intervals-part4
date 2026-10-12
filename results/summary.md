# Numerical facts behind Part IV (q = 2^k, k odd)

All numbers below are produced by the scripts in `code/`; the raw outputs are in `data/`.

## Hexad formula (Theorem 3.3), `part4_hexad_formula.py`

`s_{2m} = q^m Σ_{w ∈ T'^{-1}(1) ∩ F_{q^m}} ε(w) i^{q^{(m)}(w)}` holds for every λ1 and every pair of every class with N0 | 4 (m = 2) and N0 ∈ {3, 6} (m = 3) at q = 8, and for three pairs per class with N0 | 4 at q = 32 (m = 2). Number of points of the coset in F_{q^m}: 16 (classes 1⁶, 2 1⁴, 2² 1², 2³, 3 1³, 3²), 8 (3 2 1), 4 (4 1², 4 2, 6).

## Trace parity (Theorem 3.5, Corollary 3.6), `part4_trace_parity.py`

Traces of all rational roots, by class (f1 = dim K^σ):

| class | f1 | Tr(r) | q=8 | q=32 | q=128 |
|---|---|---|---|---|---|
| 1⁶ | 4 | 1 | – | 6 | 126 |
| 2 1⁴ | 3 | 0 | – | 60 | 1260 |
| 2² 1² | 2 | 1 | 6 | 120 | 2016 |
| 3 1³ | 2 | 1 | 9 | 165 | 2709 |
| 3 2 1 | 1 | 0 | 9 | 165 | 2709 |
| 4 1² | 1 | 0 | 12 | 240 | 4032 |
| 5 1 | 0 | 1 | 13 | 205 | 3277 |

No exception in any case.

## Classes of order 4 (Theorem 3.9, Proposition 3.12), `part4_order4_all.py`, `part4_kappa_scan.py`

* sign of s_4 = (−1)^{Tr(ρ0/τ0²)}: +1 on all 4 1² pairs, −1 on all 4 2 pairs (q = 8: all 12 pairs; q = 32: all 240 pairs, each with all 32 values of λ1; `data/part4_order4_k5.txt`), W'/q^7 = +2 and −2 on every pair.
* s_4(λ1) ≠ 0 exactly on the coset Tr(λ1 τ0) = κ0, with κ0 = ½ q^{(4)}_0(r1) = closed form Tr(z³ + ρ0²/z + (z+1)τ0²/ρ0); the trace expression (eq. 3.8) and the closed form agree with the direct Witt value on every pair.
* κ0 by class (mixed+pure pairs with κ0 = 0 / κ0 = 1): q=8: 4 1² 6/0, 4 2 6/0; q=32: 4 1² 40/80, 4 2 70/50; q=128: 4 1² 924/1092, 4 2 1106/910.
* ½ q^{(2)}_0(r5) + κ0 = 1 in the class 4 2 (both derivations); Tr(r5) = 0 in the class 4 1².
* Resolvent criterion of Lemma 3.14 (irreducible quartic cofactor ⇔ unique rational root r of Y³ + a'Y + b' with Tr(c'/(a'+r²)²) = 1, Tr(a'/r²) = 0) holds for all (y, s) at q = 8, 32, 128.

## Sign rules for the other classes (Lemma 3.16, Theorems 3.17, 3.20), `part4_sign_rules.py`

All pairs at q = 8, all λ1, and three pairs per class with N0 | 2 at q = 32 (including the classes 1⁶ and 2 1⁴, absent at q = 8): the Witt values on the coset have the predicted parity (even for m = 2, odd for m = 3), Σφ_i = 0, rational roots have φ = f1 + 1, and s_{2m}/q^m equals +4·i^p if n1 ≡ 2 (mod 4) and −4·i^p otherwise, without exception. Hypothesis (H) (q^{(3)}_0 = 0 on K ∩ F_{q³}) holds for every pair of the classes 3 1³, 3², 6 at q = 8; at q = 32 s_3 ≠ 0 for all 128 sampled triples of the classes 3² and 6 (Part III data).

## The classes 3 2 1 and 6 (Theorem 3.23), `part4_class321.py`

S'' computed from the L-polynomial of the Witt sheaf (s_1..s_4), independently of the hexad formula. Class 3 2 1: S''(λ1)·(−1)^{Tr(λ1 t0)} is constant (t0 = r4 + r5, the nonzero element of ker T' ∩ F_q) for all 9 pairs at q = 8 (constant −4q⁶) and the first 12 pairs at q = 32 (−4q⁶ for 8 pairs, +4q⁶ for 4); exactly q/2 values of λ1 have S'' = +4q⁶ in every case. Class 6: S'' = +4q⁶ for every λ1 on all 12 mixed pairs at q = 8 and the first 12 at q = 32.

## The 11 × 11 moment system, `part4_moment_system.py`

Rows m = 1..6 plus N1N2, N1N3, N2², N1²N2, N1N2N3: rank 10, det 0. Rows m = 1..6 plus N3², N2N3, N2², N1N3, N3³: rank 11, det 1433272320 = 2¹⁷·3⁷·5. Predicted mixed moments 8q²−12q−7, 5q²−11q+3, 7q²−18q+9, 3q²−7q+3, 38q²−54q−43; back-substitution returns Table 2 of Part III.

## Even k (remark after Theorem 4.1), `part4_trace_parity.py 2 4 6`, `part4_order4_all.py 2`

* Tr(r) ≡ f1 (mod 2) for every rational root at q = 4, 16, 64 (e.g. q=64: 1024 roots of class 4 1² all trace 1, 819 roots of class 5 1 all trace 0); no exception.
* q = 4: the four order-4 pairs have s_4 ∈ {0, −4q²} (class 4 1², W' = −2q^7) and {0, +4q²} (class 4 2, W' = +2q^7), on the coset Tr(λ1 τ0) = κ0 (κ0 = 0 for all four); both κ0 formulas agree with the direct Witt value.

## Moments (Section 2), `moment3_check.py`, `moment3_structure.py`, `part4_moment_system.py`

* Σ N_3 = (q−1)² + q^{−1} Σ_Θ G_H(θ): exact at q = 8, 32, 128 (103, 1951, 32383 = 2q² − 3q − 1).
* G_H(θ) ∈ {0, 2^{k+1}}, number of nonzero values (q+1)(q−2)/2 = 27, 495, 8127; dim R_H = 2 always; Arf = 0 whenever nonzero.
* Per Gold type (G_H ≠ 0 / = 0): q=8: (1,4) 27/27, (2,3) 0/9, (1,1,1,2) 0/0; q=32: 285/255, 90/243, 120/30; q=128: 4347/4095, 2436/2901, 1344/1260 — not a common quadratic polynomial per type.
* Rank of the fixed-point matrix (m = 1..6, eleven classes): 6. Predicted Σ N_m: (q−1)², 2q²−4q+1, 2q²−3q−1, 3q²−6q+1, 2q²−2q+2, 4q²−4q−1; Gauss-sum totals Σ_m: 0, q(q−1), q(q+1)(q−2), q³(2q−3), q³(q²+1), q⁴(q−1)(3q+2).
