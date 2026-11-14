"""
moment_higher_m.py -- Part IV, Section 2.3: the Gauss-sum totals Sigma_m = sum_{Theta_m} G_{H_m}(theta) for m = 4, 5, 6
and the structure that the proof of Theorem 2.2 (m = 3) does and does not provide for them.

For theta = 1/c, c in H_m = ker Tr_{q^m/q}, the F_q-line L = F_q c is totally singular for Q(e) = Tr(c^{-1}(e^5 + e^3))
(Tr(c^{-1} x^5 c^5) = Tr_{q/2}(x^5 tau_m(c^4)) = 0).  Summing over L first gives
    G_{H_m}(1/c) = q * sum_{w in K / L} (-1)^{Q(w)},   K = {w in H_m : tau_m(c^3 w) = tau_m(c w)^2}  (the orthogonal of L),
an F_2-space K/L of dimension dim K - k >= (m-3) k.  For m = 3 this is a one-dimensional space (two points), which is
the whole proof of Theorem 2.2; for m >= 4 the residual sum is a genuine Gauss sum in >= k variables.  The script
tabulates, for every theta in Theta_m, the value G_{H_m}(theta), the dimension of K/L, and the ray structure
(G along {v c : v in F_q^x}), and compares Sigma_m with the totals predicted by the class counts:
    Sigma_4 = q^3 (2q-3),  Sigma_5 = q^3 (q^2+1),  Sigma_6 = q^4 (q-1)(3q+2)        (k odd).
Usage: python moment_higher_m.py 3 4 5 6      (k, then the list of m; m = 6 at q = 8 takes a few minutes)
"""
import os
import sys
import time
from collections import Counter

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import legendre_var_check as LV  # noqa: E402
import legendre_d6_qf as QF  # noqa: E402

LV.PRIMITIVE.update({3: 0b1011, 5: 0b100101, 7: 0b10000011, 9: 0b1000010001, 15: 0b1000000000000011,
                     21: 0b1000000000000000000101, 30: 0b1000000000000000000000000000011})

PRED = {3: lambda q: q * (q + 1) * (q - 2), 4: lambda q: q ** 3 * (2 * q - 3), 5: lambda q: q ** 3 * (q * q + 1),
        6: lambda q: q ** 4 * (q - 1) * (3 * q + 2)}
PRED_N = {3: lambda q: 2 * q * q - 3 * q - 1, 4: lambda q: 3 * q * q - 6 * q + 1, 5: lambda q: 2 * q * q - 2 * q + 2,
          6: lambda q: 4 * q * q - 4 * q - 1}

OUT = open(os.path.join(HERE, "..", "data", "moment_higher_m.txt"), "a", encoding="utf-8")


def say(*a):
    s = " ".join(str(x) for x in a)
    print(s)
    OUT.write(s + "\n")
    OUT.flush()


def run(k, m):
    t0 = time.time()
    q = 1 << k
    n = m * k
    gf = LV.GF(n)
    exp, log, order, size = gf.exp, gf.log, gf.order, gf.size

    def pw(v, e):
        out = np.zeros_like(v)
        mm = v != 0
        out[mm] = exp[(e * log[v[mm]]) % order]
        return out

    def mul(u, v):
        out = np.zeros_like(u)
        mm = (u != 0) & (v != 0)
        out[mm] = exp[(log[u[mm]] + log[v[mm]]) % order]
        return out

    def tau(v):
        acc = v.copy()
        w = v
        for _ in range(m - 1):
            w = pw(w, q)
            acc = acc ^ w
        return acc

    def tr2(v):
        return LV._abs_trace_parity(gf, v.astype(np.uint32)).astype(np.int64)

    def bits(a):
        return np.array([(a >> i) & 1 for i in range(n)], dtype=np.uint8)

    allv = np.arange(size, dtype=np.int64)
    H = allv[tau(allv) == 0]
    assert len(H) == q ** (m - 1)
    Hnz = H[H != 0]
    e5e3 = np.zeros_like(H)
    e5e3[H != 0] = pw(Hnz, 5) ^ pw(Hnz, 3)
    c_all = Hnz                                              # c = 1/theta runs over H \ 0
    Fq_nz = exp[np.arange(q - 1) * (order // (q - 1))]
    # F_2-basis of H
    Mtau = np.array([bits(int(tau(np.array([1 << i], dtype=np.int64))[0])) for i in range(n)], dtype=np.uint8).T
    Hb = QF.f2_nullspace(Mtau)
    assert len(Hb) == (m - 1) * k
    Hb_el = np.array([int(sum(int(b) << i for i, b in enumerate(v))) for v in Hb], dtype=np.int64)

    vals = Counter()
    dimKL = Counter()
    ray_patterns = Counter()
    Sigma = 0
    seen = set()
    for ci in c_all:
        ci = int(ci)
        th = int(exp[(-log[ci]) % order])
        prod = mul(np.full_like(e5e3, th), e5e3)
        G = int(np.sum(1 - 2 * tr2(prod)))
        vals[G] += 1
        Sigma += G
        # dim K: K = {w in H : tau(c^3 w) + tau(c w)^2 = 0}, F_2-linear in w
        c3 = int(exp[(3 * log[ci]) % order])
        img = []
        for b in Hb_el:
            w = np.array([int(b)], dtype=np.int64)
            t1 = tau(mul(np.array([c3]), w))
            t2 = tau(mul(np.array([ci]), w))
            img.append(bits(int(t1[0] ^ pw(t2, 2)[0])))
        M = np.array(img, dtype=np.uint8).T
        dimK = len(QF.f2_nullspace(M))
        dimKL[dimK - k] += 1
    say(f"k={k}, q={q}, m={m}: |H_m| = {len(H)}, |Theta_m| = {len(c_all)}")
    say(f"   values of G_{{H_m}}(theta) -> multiplicity: {dict(sorted(vals.items()))}")
    say(f"   dim_F2 (K/L) -> multiplicity: {dict(sorted(dimKL.items()))}   (K/L is the residual space after the"
        f" totally singular line L = F_q c; for m = 3 it is {{1}})")
    pred = PRED[m](q)
    say(f"   Sigma_m = {Sigma};  predicted by the class counts: {pred};  equal: {Sigma == pred}")
    Nm = (q - 1) ** 2 - q * (m % 2 == 0) + Sigma // q ** (m - 2) if Sigma % q ** (m - 2) == 0 else None
    say(f"   sum N_m = (q-1)^2 - q[2|m] + q^(2-m) Sigma_m = {Nm}  (predicted {PRED_N[m](q)})")
    say(f"   [{time.time() - t0:.1f} s]")


def rays(k, m, max_rays=6):
    """print the multiset of G_{H_m} along a few rays {v c} and the F_q-scaling law."""
    q = 1 << k
    n = m * k
    gf = LV.GF(n)
    exp, log, order, size = gf.exp, gf.log, gf.order, gf.size

    def pw(v, e):
        out = np.zeros_like(v)
        mm = v != 0
        out[mm] = exp[(e * log[v[mm]]) % order]
        return out

    def mul(u, v):
        out = np.zeros_like(u)
        mm = (u != 0) & (v != 0)
        out[mm] = exp[(log[u[mm]] + log[v[mm]]) % order]
        return out

    def tau(v):
        acc = v.copy()
        w = v
        for _ in range(m - 1):
            w = pw(w, q)
            acc = acc ^ w
        return acc

    def tr2(v):
        return LV._abs_trace_parity(gf, v.astype(np.uint32)).astype(np.int64)

    allv = np.arange(size, dtype=np.int64)
    H = allv[tau(allv) == 0]
    Hnz = H[H != 0]
    e5e3 = np.zeros_like(H)
    e5e3[H != 0] = pw(Hnz, 5) ^ pw(Hnz, 3)
    Fq_nz = exp[np.arange(q - 1) * (order // (q - 1))]
    seen = set()
    shown = 0
    say(f"k={k}, q={q}, m={m}: G_{{H_m}} along rays {{v c : v in F_q^x}} (first {max_rays} rays)")
    for ci in Hnz:
        ci = int(ci)
        if ci in seen:
            continue
        pts = [int(exp[(log[ci] + log[v]) % order]) for v in Fq_nz]
        seen.update(pts)
        Gs = []
        for p in pts:
            th = int(exp[(-log[p]) % order])
            prod = mul(np.full_like(e5e3, th), e5e3)
            Gs.append(int(np.sum(1 - 2 * tr2(prod))))
        say(f"   ray of c={ci}: sorted G values {sorted(Gs)}  (multiset {dict(sorted(Counter(Gs).items()))})")
        shown += 1
        if shown >= max_rays:
            break


if __name__ == "__main__":
    args = [int(a) for a in sys.argv[1:]]
    k = args[0] if args else 3
    ms = args[1:] or [4, 5]
    OUT.write("\n=== moment_higher_m.py run ===\n")
    for m in ms:
        run(k, m)
        if m <= 5:
            rays(k, m)
