"""
moment3_proof.py -- Part IV, Theorem 2.2 (the third moment): pointwise check of the closed form.

For q = 2^k, tau = Tr_{q^3/q}, H = ker tau, and c in H \\ {0} with minimal polynomial c^3 + s c + N (s, N in F_q;
s != 0 when k is odd), the theorem says that the Gauss sum of e -> Tr_{q^3/2}(c^{-1}(e^5 + e^3)) on H is
    G_H(1/c) = q (1 + (-1)^{Tr_{q/2}(s^5 / N^2)})      (s != 0),
    G_H(1/c) = q                                        (s = 0, only possible for even k),
and that on every F_q^x-ray {v c} the quantity s^5/N^2 is multiplied by v^4, so each of the q + 1 rays with s != 0
contains exactly (q - 2)/2 points with G_H = 2q.  Hence
    Sigma_3 = sum_{Theta} G_H = q (q+1)(q-2)            (k odd)      -> sum N_3 = 2q^2 - 3q - 1,
    Sigma_3 = 2(q-1) q + (q-1)(q-2) q = q^2 (q-1)        (k even)     -> sum N_3 = 2q^2 - 3q + 1.
The script compares the closed form with the brute-force sum for every c, checks the ray count, and (for q <= 16)
recomputes sum N_3 by factoring the sextics.  Usage: python moment3_proof.py 3 5 7 2 4 6
"""
import os
import sys
import time
from collections import Counter

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import legendre_var_check as LV  # noqa: E402

LV.PRIMITIVE.update({3: 0b1011, 5: 0b100101, 7: 0b10000011, 9: 0b1000010001, 15: 0b1000000000000011,
                     21: 0b1000000000000000000101})

OUT = open(os.path.join(HERE, "..", "data", "moment3_proof.txt"), "a", encoding="utf-8")


def say(*a):
    s = " ".join(str(x) for x in a)
    print(s)
    OUT.write(s + "\n")
    OUT.flush()


def run(k):
    t0 = time.time()
    q = 1 << k
    n = 3 * k
    gf = LV.GF(n)
    exp, log, order, size = gf.exp, gf.log, gf.order, gf.size

    def pw(v, e):                      # array, nonzero entries
        return exp[(e * log[v]) % order]

    def mul(u, v):
        out = np.zeros_like(u)
        m = (u != 0) & (v != 0)
        out[m] = exp[(log[u[m]] + log[v[m]]) % order]
        return out

    def frq(v):
        out = np.zeros_like(v)
        m = v != 0
        out[m] = exp[(q * log[v[m]]) % order]
        return out

    def tau(v):
        v1 = frq(v)
        return v ^ v1 ^ frq(v1)

    def tr2(v):
        return LV._abs_trace_parity(gf, v.astype(np.uint32)).astype(np.int64)

    allv = np.arange(size, dtype=np.int64)
    H = allv[tau(allv) == 0]
    assert len(H) == q * q
    Hnz = H[H != 0]
    e5e3 = np.zeros_like(H)
    e5e3[H != 0] = pw(Hnz, 5) ^ pw(Hnz, 3)

    # minimal polynomial data of c in H \ 0:  s = c c^q + c c^{q^2} + c^q c^{q^2},  N = c^{1+q+q^2}
    c = Hnz
    c1 = frq(c)
    c2 = frq(c1)
    s = mul(c, c1) ^ mul(c, c2) ^ mul(c1, c2)
    N = mul(mul(c, c1), c2)
    assert np.all(frq(s) == s) and np.all(frq(N) == N) and np.all(N != 0)
    inF = np.all(frq(c) == c)
    assert not inF, "some c in H lies in F_q"
    # closed form
    phase = np.zeros_like(c)
    sn = s != 0
    phase[sn] = exp[(5 * log[s[sn]] - 2 * log[N[sn]]) % order]        # s^5 / N^2
    closed = np.where(sn, q * (1 + 1 - 2 * tr2(phase)), q)
    # brute force
    brute = np.zeros_like(c)
    for i, ci in enumerate(c):
        th = int(exp[(-log[ci]) % order])
        prod = mul(np.full_like(e5e3, th), e5e3)
        brute[i] = int(np.sum(1 - 2 * tr2(prod)))
    ok = np.array_equal(brute, closed)
    say(f"k={k}, q={q}: {len(c)} values of theta; closed form == brute force for all theta: {ok}")
    if not ok:
        bad = np.where(brute != closed)[0][:5]
        say("   first mismatches (c, brute, closed):", [(int(c[i]), int(brute[i]), int(closed[i])) for i in bad])
    say(f"   #c with s = 0: {int(np.sum(~sn))} (expected 0 for odd k, 2(q-1) for even k)")
    say(f"   values of G_H -> multiplicity: {dict(sorted(Counter(brute.tolist()).items()))}")
    Sigma = int(brute.sum())
    nz2q = int(np.sum(brute == 2 * q))
    # rays: c and v c (v in F_q^x) are on the same ray; identify rays by c^{q-1}... use c / lexicographically smallest
    Fq_nz = exp[np.arange(q - 1) * (order // (q - 1))]
    seen = {}
    ray_counts = Counter()
    ray_s0 = 0
    cset = {int(x): i for i, x in enumerate(c)}
    for i, ci in enumerate(c):
        ci = int(ci)
        if ci in seen:
            continue
        pts = [int(exp[(log[ci] + log[v]) % order]) for v in Fq_nz]
        for p in pts:
            seen[p] = True
        idx = [cset[p] for p in pts]
        if np.any(s[idx] == 0):
            assert np.all(s[idx] == 0)
            ray_s0 += 1
            ray_counts["s=0: all G_H = q" if np.all(brute[idx] == q) else "s=0: ???"] += 1
        else:
            ray_counts[int(np.sum(brute[idx] == 2 * q))] += 1
    say(f"   rays (projective points of H): {len(ray_counts) and sum(ray_counts.values())} = q+1 = {q + 1};"
        f" rays with s = 0: {ray_s0}; #points with G_H = 2q per ray -> number of rays: {dict(ray_counts)}"
        f"  (expected (q-2)/2 = {(q - 2) // 2} on every ray with s != 0)")
    if k % 2:
        pred_S, pred_N3 = q * (q + 1) * (q - 2), 2 * q * q - 3 * q - 1
    else:
        pred_S, pred_N3 = q * q * (q - 1), 2 * q * q - 3 * q + 1
    N3 = (q - 1) ** 2 + Sigma // q
    say(f"   #theta with G_H = 2q: {nz2q};  Sigma_3 = {Sigma} (predicted {pred_S}: {Sigma == pred_S});"
        f"  sum N_3 = (q-1)^2 + Sigma_3/q = {N3} (predicted {pred_N3}: {N3 == pred_N3})")
    if q <= 16:
        from moment3_check import left_side
        L = left_side(k)
        say(f"   sum N_3 by factoring all (a,b): {L}  (equal: {L == N3})")
    say(f"   [{time.time() - t0:.1f} s]")


if __name__ == "__main__":
    OUT.write("\n=== moment3_proof.py run ===\n")
    for k in [int(a) for a in sys.argv[1:]] or [3, 5]:
        run(k)
