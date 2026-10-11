# Numerical facts behind Part IV (q = 2^k, k odd)

All numbers below are produced by the scripts in `code/`; the raw outputs are in `data/`.

## Hexad formula (Theorem 3.3), `part4_hexad_formula.py`

`s_{2m} = q^m Î£_{w âˆˆ T'^{-1}(1) âˆ© F_{q^m}} Îµ(w) i^{q^{(m)}(w)}` holds for every Î»1 and every pair of every class with N0 | 4 (m = 2) and N0 âˆˆ {3, 6} (m = 3) at q = 8, and for three pairs per class with N0 | 4 at q = 32 (m = 2). Number of points of the coset in F_{q^m}: 16 (classes 1â¶, 2 1â´, 2Â² 1Â², 2Â³, 3 1Â³, 3Â²), 8 (3 2 1), 4 (4 1Â², 4 2, 6).

## Trace parity (Theorem 3.5, Corollary 3.6), `part4_trace_parity.py`

Traces of all rational roots, by class (f1 = dim K^Ïƒ):

| class | f1 | Tr(r) | q=8 | q=32 | q=128 |
|---|---|---|---|---|---|
| 1â¶ | 4 | 1 | â€“ | 6 | 126 |
| 2 1â´ | 3 | 0 | â€“ | 60 | 1260 |
| 2Â² 1Â² | 2 | 1 | 6 | 120 | 2016 |
| 3 1Â³ | 2 | 1 | 9 | 165 | 2709 |
| 3 2 1 | 1 | 0 | 9 | 165 | 2709 |
| 4 1Â² | 1 | 0 | 12 | 240 | 4032 |
| 5 1 | 0 | 1 | 13 | 205 | 3277 |

No exception in any case.

## Classes of order 4 (Theorem 3.9, Proposition 3.12), `part4_order4_all.py`, `part4_kappa_scan.py`

* sign of s_4 = (âˆ’1)^{Tr(Ï0/Ï„0Â²)}: +1 on all 4 1Â² pairs, âˆ’1 on all 4 2 pairs (q = 8: 12 pairs, all Î»1; q = 32: see `data/part4_order4_k5.txt`), W'/q^7 = +2 and âˆ’2.
* s_4(Î»1) â‰  0 exactly on the coset Tr(Î»1 Ï„0) = Îº0, with Îº0 = Â½ q^{(4)}_0(r1) = closed form Tr(zÂ³ + Ï0Â²/z + (z+1)Ï„0Â²/Ï0); the trace expression (eq. 3.8) and the closed form agree with the direct Witt value on every pair.
* Îº0 by class (mixed+pure pairs with Îº0 = 0 / Îº0 = 1): q=8: 4 1Â² 6/0, 4 2 6/0; q=32: 4 1Â² 40/80, 4 2 70/50; q=128: 4 1Â² 924/1092, 4 2 1106/910.
* Â½ q^{(2)}_0(r5) + Îº0 = 1 in the class 4 2 (both derivations); Tr(r5) = 0 in the class 4 1Â².
* Resolvent criterion of Lemma 3.14 (irreducible quartic cofactor â‡” unique rational root r of YÂ³ + a'Y + b' with Tr(c'/(a'+rÂ²)Â²) = 1, Tr(a'/rÂ²) = 0) holds for all (y, s) at q = 8, 32, 128.

## Sign rules for the other classes (Lemma 3.16, Theorems 3.17, 3.20), `part4_sign_rules.py`

All pairs at q = 8, all Î»1 (and three pairs per class with N0 | 2 at q = 32): the Witt values on the coset have the predicted parity (even for m = 2, odd for m = 3), Î£Ï†_i = 0, rational roots have Ï† = f1 + 1, and s_{2m}/q^m equals +4Â·i^p if n1 â‰¡ 2 (mod 4) and âˆ’4Â·i^p otherwise, without exception. Hypothesis (H) (q^{(3)}_0 = 0 on K âˆ© F_{qÂ³}) holds for every pair of the classes 3 1Â³, 3Â², 6 at q = 8; at q = 32 s_3 â‰  0 for all 128 sampled triples of the classes 3Â² and 6 (Part III data).

## The 11 Ã— 11 moment system, `part4_moment_system.py`

Rows m = 1..6 plus N1N2, N1N3, N2Â², N1Â²N2, N1N2N3: rank 10, det 0. Rows m = 1..6 plus N3Â², N2N3, N2Â², N1N3, N3Â³: rank 11, det 1433272320 = 2Â¹â·Â·3â·Â·5. Predicted mixed moments 8qÂ²âˆ’12qâˆ’7, 5qÂ²âˆ’11q+3, 7qÂ²âˆ’18q+9, 3qÂ²âˆ’7q+3, 38qÂ²âˆ’54qâˆ’43; back-substitution returns Table 2 of Part III.

## Even k (remark after Theorem 4.1), `part4_trace_parity.py 2 4 6`, `part4_order4_all.py 2`

* Tr(r) â‰¡ f1 (mod 2) for every rational root at q = 4, 16, 64 (e.g. q=64: 1024 roots of class 4 1Â² all trace 1, 819 roots of class 5 1 all trace 0); no exception.
* q = 4: the four order-4 pairs have s_4 âˆˆ {0, âˆ’4qÂ²} (class 4 1Â², W' = âˆ’2q^7) and {0, +4qÂ²} (class 4 2, W' = +2q^7), on the coset Tr(Î»1 Ï„0) = Îº0 (Îº0 = 0 for all four); both Îº0 formulas agree with the direct Witt value.

## Moments (Section 2), `moment3_check.py`, `moment3_structure.py`, `part4_moment_system.py`

* Î£ N_3 = (qâˆ’1)Â² + q^{âˆ’1} Î£_Î˜ G_H(Î¸): exact at q = 8, 32, 128 (103, 1951, 32383 = 2qÂ² âˆ’ 3q âˆ’ 1).
* G_H(Î¸) âˆˆ {0, 2^{k+1}}, number of nonzero values (q+1)(qâˆ’2)/2 = 27, 495, 8127; dim R_H = 2 always; Arf = 0 whenever nonzero.
* Per Gold type (G_H â‰  0 / = 0): q=8: (1,4) 27/27, (2,3) 0/9, (1,1,1,2) 0/0; q=32: 285/255, 90/243, 120/30; q=128: 4347/4095, 2436/2901, 1344/1260 â€” not a common quadratic polynomial per type.
* Rank of the fixed-point matrix (m = 1..6, eleven classes): 6. Predicted Î£ N_m: (qâˆ’1)Â², 2qÂ²âˆ’4q+1, 2qÂ²âˆ’3qâˆ’1, 3qÂ²âˆ’6q+1, 2qÂ²âˆ’2q+2, 4qÂ²âˆ’4qâˆ’1; Gauss-sum totals Î£_m: 0, q(qâˆ’1), q(q+1)(qâˆ’2), qÂ³(2qâˆ’3), qÂ³(qÂ²+1), qâ´(qâˆ’1)(3q+2).

