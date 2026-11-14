# Prime polynomials in Legendre intervals over F_{2^k}[t], IV

**Point counts on hexad modular curves, Brown invariants of Witt sheaves, and the constant X**

Ruqing Chen (GUT Geoservice Inc., Montreal) — version 2, October 2026. DOI [10.5281/zenodo.23308682](https://doi.org/10.5281/zenodo.23308682) (version 1: 10.5281/zenodo.23299284).

**New in version 3 (November 2026):** the third moment is proved. For c = 1/θ in the trace-zero plane of F_{q³} with minimal polynomial c³ + s c + N over F_q, the Gauss sum of e ↦ Tr(θ(e⁵ + e³)) on that plane is `G_H(θ) = q (1 + (−1)^{Tr(s⁵/N²)})` (Theorem 2.2); the phase s⁵/N² has weight 4 along the F_q^×-rays, so exactly (q−2)/2 points of each of the q+1 rays contribute 2q, and `Σ_λ N_3(B_λ) = 2q² − 3q − 1` for odd k (Theorem 2.3; for even k the value is 2q² − 3q + 1). The former Proposition 2.2 (numerical) and Conjecture 2.4 are replaced by these theorems; the Gold quintic does not enter. Remark 2.5 explains why the argument stops at m = 3 and confirms the predicted totals Σ_4, Σ_5, Σ_6 at q = 8 by direct summation (`code/moment_higher_m.py`). **The class counts are proved** (Theorem 3.16, Section 3.6): Lemma 3.14 parametrises the sextics with a quadratic factor by a plane (y, s) on which the quartic cofactor has resolvent Y³ + a'Y + a'/s; a rational root r is linearised by u = rs, ξ = 1/r, and the number of rational roots, the splitting of the quadratic factor and the splitting of the partition quadratic are Tr(u/(u+1)), Tr(uξ), Tr((u+1)ξ), equidistributed over the quadrants because u and u+1 are F_2-independent. The eight classes with a quadratic factor are read off cell by cell; 5 1, 3², 6 follow from the first moment, the third moment and the total. All eleven counts of Part III, Conjecture 5.3, are theorems, and **X = q⁶ + q⁵ is unconditional** (Theorem 4.1, Corollary 1.3); the higher Gauss-sum totals Σ_4, Σ_5, Σ_6 and the mixed moments follow (Corollary 2.9). Every step is checked on the (y, s)-plane at q = 8, 32, 128, 512 (`code/part4_counts_proof.py`).

**Even k, and X for all k (version 3, 13 November 2026).** The same plane proves the class counts for even k (Theorem 3.17): only Tr(1) changes, the one-root and three-root loci have q/2 and (q−4)/2 values of u, and the totals are

| 1⁶ | 2 1⁴ | 2² 1² | 2³ | 3 1³ | 3² | 3 2 1 | 4 1² | 4 2 | 5 1 | 6 |
|---|---|---|---|---|---|---|---|---|---|---|
| (q−4)(q−16)/720 | q(q−4)/48 | q(q−4)/16 | q(q−4)/48 | (q−1)(q−4)/18 | (q−1)(q+2)/18 | q(q−1)/6 | q²/8 | q²/8 | (q²−1)/5 | q(q−1)/6 |

They **do not depend on k mod 4**. The line λ3 = 0 is the curve a' = s/(s²+1) of the plane, on which the resolvent has the root r = 1/(s+1) (Proposition 3.18). This gives the pure counts for both parities: for odd k a new proof of Part III, Proposition 5.2 (4 1²: (q−2)/2, 2³: (q−2)/6, 6: (q+1)/3); for even k, 1⁶: (q−4)/60 or (q−16)/60, 2² 1²: (q−4)/4 or q/4, 5 1: 2(q+1)/5 or 2(q−1)/5 for k ≡ 2 or 0 mod 4, and 3²: (q−1)/3. So the mixed counts of Part III, Table 9, are theorems for both residues, and the whole k mod 4 dependence comes from the two points s ∈ F_4 \ F_2 missing from the pure curve, whose trace is k/2 mod 2 (Remark 3.19). The classes 3² and 6 on the pure line are separated by the F_4-structure of the kernel (σ is F_4-linear iff k is even). The sign law holds for every k with all values of W' multiplied by (−1)^{k+1} (Theorem 3.30; the classes 3 1³, 3², 5 1 now come straight from the hexad formula with m = 6). For λ5 = 0 a four-point "tetrad" analogue of the hexad formula gives W'(λ3, 0) = (−1)^{k+1}·2q⁷ or 0 including the sign (Proposition 3.32); earlier versions took this sign from the numerical data of Part II (k ≤ 7). Result: **X = (−1)^{k+1} q⁶ + q⁵ for all k** (Theorem 4.1), the formula conjectured in Part II (Conjecture 6.5). Checks: `code/part4_counts_proof.py 2 4 6 8` (every point of the plane, q = 4, 16, 64, 256), `code/part4_pure_cubic.py` (tetrad formula by direct summation, q = 2, 4, 8, 16; W'(λ3, 0) directly at q = 2, 4).

**New in version 2:** Theorem 3.28 (numbering of version 3) proves the sign laws in the classes 3 2 1 and 6, so the sign conjecture of Part III holds in all eleven classes (Theorem 1.2). X = q⁶ + q⁵ then depended only on the class counts (Theorem 1.1, proved in version 3).

Part I: DOI 10.5281/zenodo.23225896 · Part II: DOI 10.5281/zenodo.23252762 · Part III: DOI 10.5281/zenodo.23290724 ([repository](https://github.com/Ruqing1963/legendre-intervals-part3)).

## What this paper proves

Part III reduced the closed form `X = q^6 + q^5` (q = 2^k, k odd) of the second constant of the degree-6 Legendre variance to two statements about the sextic family `B_λ = λ5² u⁶ + λ3² u² + λ5 u + 1`: the polynomial class counts of its eleven factorisation types (Theorem 1.1) and the signs of the Gauss sums `s_4`, `s_6` of explicit Z/4-valued quadratic forms (Theorem 1.2). This paper proves both:

* **The hexad formula** (Theorem 3.3). For λ5 ≠ 0 and every m,
  `s_{2m}(λ) = q^m · Σ_{w ∈ T'^{-1}(1) ∩ F_{q^m}} ε(w) · i^{q^{(m)}_λ(w)}`,
  with ε = +1 at the six roots of B_λ and −1 at the ten synthematic totals. A Gauss sum over F_{q^{2m}} becomes a signed sum of at most sixteen fourth roots of unity indexed by the rational points of the hexad coset. Verified for every λ at q = 8 (m = 2, 3) and for all six classes with N0 | 4 at q = 32.
* **Trace parity** (Theorem 3.5). Every w in T'^{-1}(1) ∩ F_{q^m} has `Tr_{q^m/2}(w) ≡ mk + dim(K_λ ∩ F_{q^m}) (mod 2)` as soon as some s_m(λ1) ≠ 0. For m = 1 this is unconditional: a rational root of B_λ has absolute trace 0 in the classes 2 1⁴, 3 2 1, 4 1² and trace 1 in the classes 1⁶, 2² 1², 3 1³, 5 1 (verified for every (a, b) at q = 8, 32, 128).
* **The classes of order 4** (Theorem 3.9). `s_4(λ1) = 2q² (i^{q^{(2)}(r5)} − i^{q^{(2)}(w3)})`, with q^{(2)}(r5) = 0 in the class 4 1² and q^{(2)}(w3) = 0 in the class 4 2; hence s_4 ∈ {0, +4q²} resp. {0, −4q²}, `sgn s_4 = (−1)^{Tr(ρ0/τ0²)}`, and `W' = +2q^7` resp. `−2q^7`. Confirmed for all 12 pairs at q = 8 and all 240 pairs at q = 32, each with all λ1. The hyperplane parameter κ0 is given in closed form, `κ0 = Tr(z³ + ρ0²/z + (z+1)τ0²/ρ0)`, z = τ0/(λ5 ρ0); it is **not** constant on the class (130 of the 240 pairs at q = 32 have κ0 = 1), which corrects a conjecture of the first draft.
* **The other classes** (Theorems 3.22, 3.25). On the sixteen-point coset the half-Witt values are affine-linear in the duads, so the hexad sum is ±4 according to the number n₁ of roots with φ = 1 mod 4, and trace parity fixes φ at rational points. This proves `s_4 = +4q²` for all λ1 in the classes 1⁶, 2³, the half-half law in 2 1⁴, 2² 1², and `S'' = +4q⁶` for all λ1 in 3 1³, 3².
* **The classes 3 2 1 and 6** (Theorem 3.28). Applying the hexad formula with m = 6 expresses S'' = s_12 itself as a sixteen-point sum. On the coset, ε = (−1)^Q for a quadratic form Q polarising the intersection form. The λ1-free part of the Witt form is a Frobenius-invariant linear form, hence dual to a Frobenius-fixed duad. Result: `S''(λ1) = ±4q⁶ (−1)^{Tr(λ1 (r4+r5))}` in the class 3 2 1, so W' = 0, and `S'' = +4q⁶` for every λ1 in the class 6, without the hypothesis (H). **The sign conjecture of Part III (Theorem 1.2) is proved in all eleven classes.** Checked against independent L-polynomial data at q = 8 (all pairs) and q = 32 (12 pairs per class): `code/part4_class321.py`.
* **Moments** (Propositions 2.1, 2.6, 2.8; Theorems 2.2, 2.3). For every m, `Σ N_m = (q−1)² − q[2|m] + q^{2−m} Σ_{Θ_m} G_{H_m}(θ)` with G_{H_m} the Gauss sums of Gold-type forms on the trace-zero hyperplane; m = 1, 2 recover Part III. For m = 3 the Gauss sum is evaluated in closed form, `G_H(1/c) = q(1 + (−1)^{Tr(s⁵/N²)})` with c³ + s c + N the minimal polynomial of c over F_q, and the ray count gives `Σ N_3 = 2q² − 3q − 1` for odd k (2q² − 3q + 1 for even k): a theorem for every q. Checked pointwise at q = 4, 8, 16, 32, 64, 128 (`code/moment3_proof.py`). The rows m = 4, 5, 6 of the moment table follow from the class counts (Corollary 2.9).
* **The 11 × 11 moment system** (Proposition 2.8): the six pure moments plus the mixed moments Σ N3², Σ N2N3, Σ N2², Σ N1N3, Σ N3³ form a nonsingular system (determinant 2¹⁷·3⁷·5) whose unique solution is the table; the five mixed moments N1N2, N1N3, N2², N1²N2, N1N2N3 give rank 10 only.
* **The class counts** (Theorem 3.16): the eight classes with a quadratic factor by the point count on the resolvent locus of Lemma 3.14, the classes 5 1, 3², 6 from the first moment, the third moment and the total. Part III, Conjecture 5.3 is a theorem.

* **Even k and the terms λ5 = 0** (Theorems 3.17, 3.30, Propositions 3.18, 3.32): see the paragraph "Even k, and X for all k" above. The trace parity law becomes `Tr(r) ≡ f_1 (mod 2)` (verified at q = 4, 16, 64) and every class value of W' changes sign.

All inputs of X are now theorems, so **X = (−1)^{k+1} q⁶ + q⁵ holds for all k** (Theorem 4.1). Not proved here: a direct evaluation of the Gauss-sum totals Σ_m for m = 4, 5, 6 (they follow from the class counts); the identity (H), which only concerns the value of s_6 in the class 6; the constant A_6 for even k and the constant Y (so the variance formula of Corollary 1.3 is stated for odd k, with m_k = Y/q³ unevaluated); and any link between the k mod 4 branching and the fifth roots of unity on the Gold-quintic side (Section 5.2 records what is and is not proved).

## Layout

```
paper/    legendre_IV_hexad_moments_brown_X.tex / .pdf (36 pages), make_figures4.py, figures/
code/     self-contained Python (NumPy, SymPy, Numba): field tables, hexad search, Witt sums, all Part IV scripts
data/     outputs of every script quoted in the paper (q = 8, 32, 128)
figures/  the figures of the paper (vector PDF)
results/  summary.md: the numerical facts in one page
```

## Reproduce

```
cd code
python moment3_proof.py 3 5 7 2 4 6       # Theorems 2.2, 2.3: closed form of G_H, rays       -> data/moment3_proof.txt
python moment_higher_m.py 3 4 5 6         # Sigma_m, G_{H_m} values, dim K/L, m = 4,5,6, q=8 -> data/moment_higher_m.txt
python moment3_check.py 3 5 7             # third moment identity, G_H distribution         -> data/moment3_check.txt
python moment3_structure.py 3 5 7         # radicals, Arf invariants, Gold types            -> data/moment3_structure.txt
python part4_moment_system.py             # 6x11 and 11x11 moment systems, totals, fits     -> data/part4_moment_system.txt
python part4_hexad_formula.py 3           # hexad formula, all classes, q=8, m=2,3           -> data/part4_hexad_formula_k3.txt
python part4_hexad_formula.py 5 2 3       # hexad formula, q=32, m=2, 3 pairs per class      -> data/part4_hexad_formula_k5.txt
python part4_trace_parity.py 3 5 7        # traces of rational roots by class               -> data/part4_trace_parity.txt
python part4_sign_rules.py 3              # n1 sign rules and hypothesis (H), q=8           -> data/part4_sign_rules_k3.txt
python part4_sign_rules.py 5 3            # same at q=32, 3 pairs per class                 -> data/part4_sign_rules_k5.txt
python part4_kappa_scan.py 3 5 7          # kappa_0 closed form, resolvent criterion         -> data/part4_kappa_scan.txt
python part4_counts_proof.py 3 5 7 9      # Theorem 3.16: the class counts, cell by cell     -> data/part4_counts_proof.txt
python part4_counts_proof.py 2 4 6 8      # Theorem 3.17, Prop. 3.18: even k, pure curve     -> data/part4_counts_proof.txt (appended)
python part4_pure_cubic.py 1 2 3 4        # Prop. 3.32: tetrad formula, W'(lambda_3, 0)      -> data/part4_pure_cubic.txt
python part4_order4_all.py 3              # all order-4 pairs: Witt values, kappa_0, s_4     -> data/part4_order4_k3.txt
python part4_order4_all.py 5              # same at q=32 (all 240 pairs, about two hours)    -> data/part4_order4_k5.txt
python part4_relative_reduction.py 3      # coset form of the hexad formula                  -> data/part4_relative_reduction_k3.txt
cd ../paper && python make_figures4.py && pdflatex legendre_IV_hexad_moments_brown_X.tex
```

Field sizes: q = 8, 32, 128 use F_{2^{3k}} (moments), F_{2^{4k}} (s_4), F_{2^{6k}} (s_6, q = 8 only) and F_{2^{12k}} (hexad search, q ≤ 32).

## Zenodo record

Reserved DOI: [10.5281/zenodo.23308682](https://doi.org/10.5281/zenodo.23308682) (version 2; version 1: 10.5281/zenodo.23299284). Resource type: publication / preprint. Version v2, 2026-10-12, English, CC BY 4.0.

**Title.** Prime polynomials in Legendre intervals over F_{2^k}[t], IV: point counts on hexad modular curves, Brown invariants of Witt sheaves, and the constant X

**Creator.** Chen, Ruqing (GUT Geoservice Inc., Montreal).

**Description.**

> Let q = 2^k with k odd. Part III (DOI 10.5281/zenodo.23290724) reduced the closed form X = q^6 + q^5 of the second constant of the degree-6 Legendre variance to two statements about the sextic family B_λ = λ5² u⁶ + λ3² u² + λ5 u + 1, whose roots form the hexad of ker L'_λ: the polynomial frequencies of its eleven factorisation types over F_q (Theorem 1.1, Conjecture 5.3 of Part III), and the signs of the Gauss sums s_4, s_6, s_12 of explicit Z/4-valued quadratic forms (Theorem 1.2, a conjecture of Part III). Both are proved here.
>
> This paper proves three structural theorems about these Gauss sums. The hexad formula expresses s_{2m}, a sum over F_{q^{2m}}, as q^m times a signed sum of at most sixteen fourth roots of unity indexed by the F_{q^m}-rational points of the coset T'^{-1}(1) (the six roots and the ten synthematic totals). The trace parity theorem determines the absolute trace of every such point, Tr(w) ≡ mk + dim(K_λ ∩ F_{q^m}) (mod 2); for m = 1 it fixes the trace of every rational root of B_λ by Frobenius class. A sixteen-point lemma on the configuration of roots and totals then determines the signs in all eleven classes. s_4 = +4q² throughout the classes 1^6 and 2^3. A half-half law in λ1 holds in the classes 2 1^4 and 2^2 1^2. In the classes 4 1^2 and 4 2, s_4 ∈ {0, +4q²}, respectively {0, −4q²}, nonzero exactly on an explicit affine hyperplane, hence W' = +2q^7 and −2q^7 as conjectured since Part II. S'' = +4q^6 throughout the classes 3 1^3, 3^2 and 6. In the class 3 2 1, S''(λ1) = ±4q^6 (−1)^{Tr(λ1(r4+r5))}, hence W' = 0. The last two classes come from applying the hexad formula directly over F_{q^6}; version 2 adds this step, which removes the Witt identity (H) previously needed for the class 6. The sign conjecture of Part III is therefore a theorem. The hyperplane parameter κ0 of the order-4 classes is given in closed form; contrary to what the data at q = 8 suggested, it is not constant on the class.
>
> For the factorisation counts, every moment Σ_λ N_m(B_λ) is expressed through Gauss sums of Gold-type forms on the trace-zero hyperplane of F_{q^m}; the identity is proved for all m, and the third moment is evaluated in closed form: G_H(1/c) = q(1 + (−1)^{Tr(s⁵/N²)}) for c with minimal polynomial c³ + sc + N, giving Σ N_3 = 2q² − 3q − 1 for odd k. The class counts themselves are then proved: the sextics with a quadratic factor are parametrised by a plane on which the quartic cofactor has resolvent Y³ + a'Y + a'/s; a rational root r is linearised by u = rs, ξ = 1/r, and the number of rational roots, the splitting of the quadratic factor and the splitting of the partition quadratic are the traces of u/(u+1), uξ and (u+1)ξ, equidistributed because u and u+1 are F_2-independent. The eight classes with a quadratic factor follow cell by cell, the other three from the first moment, the third moment and the total. All eleven class counts of Part III, Conjecture 5.3, are theorems, and X = q^6 + q^5 holds for all odd k. The same plane gives the counts for even k: the totals depend only on the parity of k, and the dependence on k mod 4 observed in Part III lies entirely in the division between λ3 = 0 and λ3 ≠ 0, where it comes from the two points s ∈ F_4 \ F_2 missing from the pure curve a' = s/(s²+1). The sign law is proved for every k, the terms with λ5 = 0 are evaluated by a four-point analogue of the hexad formula, and X = (−1)^{k+1} q^6 + q^5 for all k, as conjectured in Part II. An explicit nonsingular 11 × 11 system of pure and mixed moments, whose unique solution is the table of class counts, is also written down.
>
> All numerical statements are reproducible with the self-contained Python code in the companion repository (https://github.com/Ruqing1963/legendre-intervals-part4), with the raw outputs for q = 8, 32, 128 included. 24 pages, 3 figures, 3 tables.
>
> Related records: Part I, DOI 10.5281/zenodo.23225896; Part II, DOI 10.5281/zenodo.23252762; Part III, DOI 10.5281/zenodo.23290724; Part IV version 1, DOI 10.5281/zenodo.23299284.

**Keywords.** Legendre intervals; prime polynomial distribution; function fields; even characteristic; Witt vectors; Z/4-valued quadratic forms; Brown invariant; Gauss sums; Arf invariant; hexad; synthematic totals; Sp4(F2); Gold-type functions; Kloosterman sums; quartic cofactor; moment method; variance of primes; finite fields.

**Subjects (MSC 2020).** 11T55 arithmetic theory of polynomial rings over finite fields; 11T24 other character sums and Gauss sums; 11T06 polynomials over finite fields; 11E81 algebraic theory of quadratic forms, Witt groups and rings; 11G20 curves over finite and local fields; 11L05 Gauss and Kloosterman sums; 14G15 finite ground fields in algebraic geometry.

**Related identifiers.** This repository (is supplemented by); Parts I–III above (continues).

## License

Text and figures: CC BY 4.0. Code: MIT. See `LICENSE`.
