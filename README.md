# Prime polynomials in Legendre intervals over F_{2^k}[t], IV

**Point counts on hexad modular curves, Brown invariants of Witt sheaves, and the constant X**

Ruqing Chen (GUT Geoservice Inc., Montreal) — working draft, October 2026.

Part I: DOI 10.5281/zenodo.23225896 · Part II: DOI 10.5281/zenodo.23252762 · Part III: DOI 10.5281/zenodo.23290724 ([repository](https://github.com/Ruqing1963/legendre-intervals-part3)).

## What this paper proves

Part III reduced the closed form `X = q^6 + q^5` (q = 2^k, k odd) of the second constant of the degree-6 Legendre variance to two statements about the sextic family `B_λ = λ5² u⁶ + λ3² u² + λ5 u + 1`: the polynomial class counts of its eleven factorisation types (Target Theorem 1.1) and the signs of the Gauss sums `s_4`, `s_6` of explicit Z/4-valued quadratic forms (Target Theorem 1.2). This paper proves:

* **The hexad formula** (Theorem 3.3). For λ5 ≠ 0 and every m,
  `s_{2m}(λ) = q^m · Σ_{w ∈ T'^{-1}(1) ∩ F_{q^m}} ε(w) · i^{q^{(m)}_λ(w)}`,
  with ε = +1 at the six roots of B_λ and −1 at the ten synthematic totals. A Gauss sum over F_{q^{2m}} becomes a signed sum of at most sixteen fourth roots of unity indexed by the rational points of the hexad coset. Verified for every λ at q = 8 (m = 2, 3) and for all six classes with N0 | 4 at q = 32.
* **Trace parity** (Theorem 3.4). Every w in T'^{-1}(1) ∩ F_{q^m} has `Tr_{q^m/2}(w) ≡ mk + dim(K_λ ∩ F_{q^m}) (mod 2)` as soon as some s_m(λ1) ≠ 0. For m = 1 this is unconditional: a rational root of B_λ has absolute trace 0 in the classes 2 1⁴, 3 2 1, 4 1² and trace 1 in the classes 1⁶, 2² 1², 3 1³, 5 1 (verified for every (a, b) at q = 8, 32, 128).
* **The classes of order 4** (Theorem 3.6). `s_4(λ1) = 2q² (i^{q^{(2)}(r5)} − i^{q^{(2)}(w3)})`, with q^{(2)}(r5) = 0 in the class 4 1² and q^{(2)}(w3) = 0 in the class 4 2; hence s_4 ∈ {0, +4q²} resp. {0, −4q²}, `sgn s_4 = (−1)^{Tr(ρ0/τ0²)}`, and `W' = +2q^7` resp. `−2q^7`. This settles the order-4 part of Target Theorem 1.2. The hyperplane parameter κ0 is given in closed form, `κ0 = Tr(z³ + ρ0²/z + (z+1)τ0²/ρ0)`, z = τ0/(λ5 ρ0); it is **not** constant on the class (130 of the 240 pairs at q = 32 have κ0 = 1), which corrects a conjecture of the first draft.
* **Moments** (Propositions 2.1, 2.6, 2.8). For every m, `Σ N_m = (q−1)² − q[2|m] + q^{2−m} Σ_{Θ_m} G_{H_m}(θ)` with G_{H_m} the Gauss sums of Gold-type forms on the trace-zero hyperplane; m = 1, 2 recover Part III; m = 3 gives 2q² − 3q − 1 (verified q ≤ 128); the elimination system of the six moments (rank 6 in the eleven unknowns) and the Gauss-sum totals predicted by the class counts are tabulated.

Even k: the trace parity law becomes `Tr(r) ≡ f_1 (mod 2)` (verified at q = 4, 16, 64) and the two order-4 signs are exchanged (`W' = −2q^7` for 4 1², `+2q^7` for 4 2; verified at q = 4), in line with the sign (−1)^{k+1} of Part II.

Not proved here (stated with precise status in the paper): the class counts; the signs of s_4 in the classes 1⁶, 2³, 2 1⁴, 2² 1² and of S'' in the classes 3 1³, 3², 3 2 1, 6 — all reduced by the hexad formula to sixteen explicit Witt values on the hexad.

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

Field sizes: q = 8, 32, 128 use F_{2^{3k}} (moments), F_{2^{4k}} (s_4), F_{2^{6k}} (s_6, q = 8 only) and F_{2^{12k}} (hexad search, q ≤ 32).

## License

Text and figures: CC BY 4.0. Code: MIT. See `LICENSE`.
