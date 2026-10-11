# Prime polynomials in Legendre intervals over F_{2^k}[t], IV

**Point counts on hexad modular curves, Brown invariants of Witt sheaves, and the constant X**

Ruqing Chen (GUT Geoservice Inc., Montreal) â€” working draft, October 2026.

Part I: DOI 10.5281/zenodo.23225896 Â· Part II: DOI 10.5281/zenodo.23252762 Â· Part III: DOI 10.5281/zenodo.23290724 ([repository](https://github.com/Ruqing1963/legendre-intervals-part3)).

## What this paper proves

Part III reduced the closed form `X = q^6 + q^5` (q = 2^k, k odd) of the second constant of the degree-6 Legendre variance to two statements about the sextic family `B_Î» = Î»5Â² uâ¶ + Î»3Â² uÂ² + Î»5 u + 1`: the polynomial class counts of its eleven factorisation types (Target Theorem 1.1) and the signs of the Gauss sums `s_4`, `s_6` of explicit Z/4-valued quadratic forms (Target Theorem 1.2). This paper proves:

* **The hexad formula** (Theorem 3.3). For Î»5 â‰  0 and every m,
  `s_{2m}(Î») = q^m Â· Î£_{w âˆˆ T'^{-1}(1) âˆ© F_{q^m}} Îµ(w) Â· i^{q^{(m)}_Î»(w)}`,
  with Îµ = +1 at the six roots of B_Î» and âˆ’1 at the ten synthematic totals. A Gauss sum over F_{q^{2m}} becomes a signed sum of at most sixteen fourth roots of unity indexed by the rational points of the hexad coset. Verified for every Î» at q = 8 (m = 2, 3) and for all six classes with N0 | 4 at q = 32.
* **Trace parity** (Theorem 3.5). Every w in T'^{-1}(1) âˆ© F_{q^m} has `Tr_{q^m/2}(w) â‰¡ mk + dim(K_Î» âˆ© F_{q^m}) (mod 2)` as soon as some s_m(Î»1) â‰  0. For m = 1 this is unconditional: a rational root of B_Î» has absolute trace 0 in the classes 2 1â´, 3 2 1, 4 1Â² and trace 1 in the classes 1â¶, 2Â² 1Â², 3 1Â³, 5 1 (verified for every (a, b) at q = 8, 32, 128).
* **The classes of order 4** (Theorem 3.9). `s_4(Î»1) = 2qÂ² (i^{q^{(2)}(r5)} âˆ’ i^{q^{(2)}(w3)})`, with q^{(2)}(r5) = 0 in the class 4 1Â² and q^{(2)}(w3) = 0 in the class 4 2; hence s_4 âˆˆ {0, +4qÂ²} resp. {0, âˆ’4qÂ²}, `sgn s_4 = (âˆ’1)^{Tr(Ï0/Ï„0Â²)}`, and `W' = +2q^7` resp. `âˆ’2q^7`. This settles the order-4 part of Target Theorem 1.2. The hyperplane parameter Îº0 is given in closed form, `Îº0 = Tr(zÂ³ + Ï0Â²/z + (z+1)Ï„0Â²/Ï0)`, z = Ï„0/(Î»5 Ï0); it is **not** constant on the class (130 of the 240 pairs at q = 32 have Îº0 = 1), which corrects a conjecture of the first draft.
* **Moments** (Propositions 2.1, 2.6, 2.8). For every m, `Î£ N_m = (qâˆ’1)Â² âˆ’ q[2|m] + q^{2âˆ’m} Î£_{Î˜_m} G_{H_m}(Î¸)` with G_{H_m} the Gauss sums of Gold-type forms on the trace-zero hyperplane; m = 1, 2 recover Part III; m = 3 gives 2qÂ² âˆ’ 3q âˆ’ 1 (verified q â‰¤ 128); the elimination system of the six moments (rank 6 in the eleven unknowns) and the Gauss-sum totals predicted by the class counts are tabulated.

Even k: the trace parity law becomes `Tr(r) â‰¡ f_1 (mod 2)` (verified at q = 4, 16, 64) and the two order-4 signs are exchanged (`W' = âˆ’2q^7` for 4 1Â², `+2q^7` for 4 2; verified at q = 4), in line with the sign (âˆ’1)^{k+1} of Part II.

* **The other classes** (Theorems 3.17, 3.20). On the sixteen-point coset the half-Witt values are affine-linear in the duads, so the hexad sum is Â±4 according to the number nâ‚ of roots with Ï† = 1 mod 4, and trace parity fixes Ï† at rational points. This proves `s_4 = +4qÂ²` for all Î»1 in the classes 1â¶, 2Â³, the half-half law in 2 1â´, 2Â² 1Â², and `S'' = +4qâ¶` for all Î»1 in 3 1Â³, 3Â². For the class 6 the same follows from one Witt identity (H) verified at q = 8, 32. Nine of the eleven sign statements of Target Theorem 1.2 are thus theorems, one is conditional, and only the class 3 2 1 is open.
* **The 11 Ã— 11 moment system** (Proposition 2.8): the six pure moments plus the mixed moments Î£ N3Â², Î£ N2N3, Î£ N2Â², Î£ N1N3, Î£ N3Â³ form a nonsingular system (determinant 2Â¹â·Â·3â·Â·5) whose unique solution is the conjectured table; the five mixed moments the referee proposed (N1N2, N1N3, N2Â², N1Â²N2, N1N2N3) give rank 10 only.

Not proved here (stated with precise status in the paper): the class counts (equivalently the eleven point counts of the moment system); the sign in the class 3 2 1; the identity (H) for the class 6.

## Layout

```
paper/    legendre_IV_hexad_moments_brown_X.tex / .pdf (20 pages), make_figures4.py, figures/
code/     self-contained Python (NumPy, SymPy, Numba): field tables, hexad search, Witt sums, all Part IV scripts
data/     outputs of every script quoted in the paper (q = 8, 32, 128)
figures/  the two figures of the paper (vector PDF)
results/  summary.md: the numerical facts in one page
```

## Reproduce

```
cd code
python moment3_check.py 3 5 7             # third moment identity, G_H distribution         -> data/moment3_check.txt
python moment3_structure.py 3 5 7         # radicals, Arf invariants, Gold types            -> data/moment3_structure.txt
python part4_moment_system.py             # elimination system, Gauss-sum totals, fits      -> data/part4_moment_system.txt
python part4_hexad_formula.py 3           # hexad formula, all classes, q=8, m=2,3           -> data/part4_hexad_formula_k3.txt
python part4_hexad_formula.py 5 2 3       # hexad formula, q=32, m=2, 3 pairs per class      -> data/part4_hexad_formula_k5.txt
python part4_trace_parity.py 3 5 7        # traces of rational roots by class               -> data/part4_trace_parity.txt
python part4_kappa_scan.py 3 5 7          # kappa_0 closed form, resolvent criterion         -> data/part4_kappa_scan.txt
python part4_order4_all.py 3              # all order-4 pairs: Witt values, kappa_0, s_4     -> data/part4_order4_k3.txt
python part4_order4_all.py 5              # same at q=32 (all 240 pairs, about an hour)      -> data/part4_order4_k5.txt
python part4_relative_reduction.py 3      # coset form of the hexad formula                  -> data/part4_relative_reduction_k3.txt
cd ../paper && python make_figures4.py && pdflatex legendre_IV_hexad_moments_brown_X.tex
```

Field sizes: q = 8, 32, 128 use F_{2^{3k}} (moments), F_{2^{4k}} (s_4), F_{2^{6k}} (s_6, q = 8 only) and F_{2^{12k}} (hexad search, q â‰¤ 32).

## License

Text and figures: CC BY 4.0. Code: MIT. See `LICENSE`.

