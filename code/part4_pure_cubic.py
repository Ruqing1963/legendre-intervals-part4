r"""part4_pure_cubic.py -- the terms lambda_5 = 0 of X: the tetrad formula (Proposition 3.27) and W'(lambda_3, 0).

For lambda_5 = 0, lambda_3 != 0 the polar operator of q_lambda(x) = T(x) + 2 Tr(lambda_3 x^3 + lambda_1 x) is
T'(x) = lambda_3 x^2 + x + (lambda_3 x)^{1/2}; its kernel K is {0} and the roots of lambda_3^2 x^3 + x + lambda_3, and
    T'^{-1}(1) = { w_0 = 1/lambda_3 }  u  { roots of lambda_3 w^3 + w^2 + 1 }   (= w_0 + K).
Tetrad formula:   s_{2m}(lambda) = q^m * sum_{w in T'^{-1}(1), w in F_{q^m}} eps(w) i^{q^(m)(w)},   eps(w_0) = +1, eps = -1 otherwise.
Consequence (m = 6):  S'' = s_12 = (-1)^{Tr(1/lambda_3)} q^6 (2 - 4 [psi = 0 on K]),  and
    W'(lambda_3, 0) = sum_{lambda_1} S'' = (-1)^{k+1} 2 q^7  if z^3 + z + lambda_3^2 has 3 or 0 roots in F_q,  0 if it has exactly one.
The script checks the tetrad formula by direct summation over F_{q^{2m}} for every (lambda_3, lambda_1), for
    q = 2: m = 1..6,   q = 4: m = 1, 2, 3, 6,   q = 8: m = 1, 2, 3, 4,   q = 16: m = 1, 2, 3,
and, where m = 6 is reachable (q = 2, 4), the values of S'' and W' themselves.
Usage: python part4_pure_cubic.py [k ...]      (default 1 2 3 4)"""
import sys, os, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = open(os.path.join(HERE, "..", "data", "part4_pure_cubic.txt"), "a", encoding="utf8")
def say(s=""):
    print(s, flush=True); OUT.write(s + "\n"); OUT.flush()

PRIM = {2: 0b111, 4: 0b10011, 6: 0b1000011, 8: 0b100011101, 10: (1 << 10) | 0b1001, 12: (1 << 12) | 0b1010011,
        16: (1 << 16) | 0b101101, 18: (1 << 18) | (1 << 7) | 1, 24: (1 << 24) | 0b10000111}
MS = {1: (1, 2, 3, 4, 5, 6), 2: (1, 2, 3, 6), 3: (1, 2, 3, 4), 4: (1, 2, 3)}

class BigGF:
    """GF(2^n) with exp/log tables; x is the generator (the modulus is checked to be primitive)."""
    def __init__(self, n):
        self.n = n; self.size = 1 << n; self.N = self.size - 1; mod = PRIM[n]
        exp = np.zeros(self.N, np.int64); x = 1
        for i in range(self.N):
            exp[i] = x; x <<= 1
            if x >> n: x ^= mod
        assert x == 1 and len(set(exp[:min(self.N, 1 << 16)].tolist())) == min(self.N, 1 << 16), "modulus not primitive"
        log = np.full(self.size, -1, np.int64); log[exp] = np.arange(self.N)
        assert (log[1:] >= 0).all(), "modulus not primitive"
        self.exp, self.log = exp, log
        # absolute trace and half trace of every element, indexed by the logarithm
        L = np.arange(self.N, dtype=np.int64)
        acc = np.zeros(self.N, np.int64); half = None; e = L.copy()
        for i in range(n):
            acc ^= exp[e]; e = (2 * e) % self.N
            if i == n // 2 - 1: half = acc.copy()
        assert set(np.unique(acc).tolist()) <= {0, 1}
        self.TRL = acc; self.HTL = half          # TRL[L] = Tr(g^L);  HTL[L] = sum_{i < n/2} (g^L)^{2^i}
    def mul(self, a, b):
        if a == 0 or b == 0: return 0
        return int(self.exp[(self.log[a] + self.log[b]) % self.N])
    def inv(self, a): return int(self.exp[(-self.log[a]) % self.N])
    def pw(self, a, e):
        if a == 0: return 0
        return int(self.exp[(self.log[a] * e) % self.N])
    def witt_all(self):
        """T_n(g^L) mod 4 for all L:  Tr + 2 e_2,  e_2 = sum_{0<d<n/2} Tr(x^{1+2^d}) + Tr_{2^{n/2}/2}(x^{1+2^{n/2}})."""
        n, N = self.n, self.N; L = np.arange(N, dtype=np.int64)
        e2 = np.zeros(N, np.int64)
        for d in range(1, n // 2):
            e2 ^= self.TRL[(L * (1 + (1 << d))) % N]
        h = self.HTL[(L * (1 + (1 << (n // 2)))) % N]
        assert set(np.unique(h).tolist()) <= {0, 1}
        e2 ^= h
        return (self.TRL + 2 * e2) % 4

def gauss(Qv):
    c = np.bincount(Qv % 4, minlength=4)
    return complex(1 + c[0] - c[2], c[1] - c[3])       # the term x = 0 contributes 1

def run(k):
    q = 1 << k
    say("\n=== k = %d, q = %d ===" % (k, q))
    for m in MS[k]:
        t0 = time.time(); n = 2 * m * k; G = BigGF(n); N = G.N; L = np.arange(N, dtype=np.int64)
        W = G.witt_all()
        stepq = N // (q - 1) if q > 2 else 0
        Fq = [0] + [int(G.exp[(j * stepq) % N]) for j in range(q - 1)] if q > 2 else [0, 1]
        M = (1 << (m * k)) - 1; stepm = N // M
        sub = [0] + [int(G.exp[j * stepm]) for j in range(M)]          # F_{q^m}
        def witt_sub(w):          # T_{mk}(w) for w in F_{q^m}
            s0 = s1 = 0; c = w
            for _ in range(m * k):
                s1 ^= G.mul(s0, c); s0 ^= c; c = G.mul(c, c)
            assert s0 in (0, 1) and s1 in (0, 1)
            return s0 + 2 * s1
        def tr_sub(w):            # Tr_{q^m/2}(w)
            t = 0; c = w
            for _ in range(m * k):
                t ^= c; c = G.mul(c, c)
            assert t in (0, 1); return t
        def trq(c):               # Tr_{q/2} of an element of F_q
            t = 0
            for _ in range(k):
                t ^= c; c = G.mul(c, c)
            assert t in (0, 1); return t
        nok = ntot = 0; rows = []
        for l3 in Fq[1:]:
            l3sq = G.mul(l3, l3); w0 = G.inv(l3)
            # coset T'^{-1}(1) in F_{q^m}: l3^2 w^4 + w^2 + l3 w + 1 = 0
            coset = [w for w in sub if G.mul(l3sq, G.pw(w, 4)) ^ G.mul(w, w) ^ G.mul(l3, w) ^ 1 == 0]
            assert w0 in coset
            for w in coset:
                if w != w0: assert G.mul(l3, G.pw(w, 3)) ^ G.mul(w, w) ^ 1 == 0
            nz = sum(1 for z in Fq if G.pw(z, 3) ^ z ^ l3sq == 0)      # rational roots of z^3 + z + l3^2
            base = (W + 2 * G.TRL[(3 * L + G.log[l3]) % N]) % 4
            Ssum = 0; vals = []
            for l1 in Fq:
                Qv = base if l1 == 0 else (base + 2 * G.TRL[(L + G.log[l1]) % N]) % 4
                lhs = gauss(Qv)
                rhs = 0
                for w in coset:
                    qw = (witt_sub(w) + 2 * tr_sub(G.mul(l3, G.pw(w, 3)) ^ G.mul(l1, w))) % 4
                    rhs += (1 if w == w0 else -1) * (1j ** qw)
                rhs *= q ** m
                ntot += 1; nok += abs(lhs - rhs) < 1e-6
                Ssum += lhs; vals.append(lhs)
            rows.append((l3, nz, len(coset), trq(G.inv(l3)), vals, Ssum))
        say("m = %d (sum over F_{2^%d}): tetrad formula holds for %d of %d pairs (lambda_3, lambda_1)   [%.0f s]" % (m, n, nok, ntot, time.time() - t0))
        if m == 6:
            okW = okS = True; hist = {}
            for (l3, nz, nc, tr1, vals, Ssum) in rows:
                pred = 0 if nz == 1 else (-1) ** (k + 1) * 2 * q ** 7
                okW &= abs(Ssum - pred) < 1e-6
                okS &= all(abs(abs(v.real) - 2 * q ** 6) < 1e-6 and abs(v.imag) < 1e-6 for v in vals)
                okS &= (nz == 1) == (tr1 == (k + 1) % 2)
                if nz != 1: okS &= all(abs(v - (-1) ** (k + 1) * 2 * q ** 6) < 1e-6 for v in vals)
                else: okS &= sum(1 for v in vals if v.real > 0) == q // 2
                hist[nz] = hist.get(nz, 0) + 1
            tot = sum(r[5] for r in rows)
            say("  S'' = s_12 = +-2q^6 for every (lambda_3, lambda_1); constant (-1)^(k+1) 2q^6 in lambda_1 when the cubic has 0 or 3 roots, "
                "half/half when it has one root (Tr(1/lambda_3) = k+1): %s" % okS)
            say("  W'(lambda_3, 0) = (-1)^(k+1) 2q^7 [0 or 3 roots], 0 [one root]: %s;  cubics by number of roots: %s" % (okW, dict(sorted(hist.items()))))
            pr = q ** 8 if k % 2 else -(q - 2) * q ** 7
            say("  sum_{lambda_3 != 0} W'(lambda_3, 0) = %d = %s: %s" % (round(tot.real), "q^8" if k % 2 else "-(q-2) q^7", abs(tot - pr) < 1e-6))

def counts():
    """number of c in F_q^x with z^3 + z + c of each root number: (q-2)/6, (q+1)/3, (q-2)/2 for odd k; (q-4)/6, (q-1)/3, q/2 for even k"""
    say("\n=== roots of z^3 + z + c, c != 0 ===")
    from part4_counts_proof import GF
    for k in (2, 3, 4, 5, 6, 7, 8, 9, 10):
        F = GF(k); q = F.q; h = {0: 0, 1: 0, 3: 0}
        img = {}
        for z in range(q):
            c = F.mul(z, F.mul(z, z)) ^ z; img[c] = img.get(c, 0) + 1
        for c in range(1, q): h[img.get(c, 0)] += 1
        pred = {3: (q - 2) // 6, 0: (q + 1) // 3, 1: (q - 2) // 2} if k % 2 else {3: (q - 4) // 6, 0: (q - 1) // 3, 1: q // 2}
        say("  q = %4d: three roots %d, none %d, one %d; predicted %d, %d, %d: %s; nonzero W' for %d = %s values"
            % (q, h[3], h[0], h[1], pred[3], pred[0], pred[1], h == pred, h[3] + h[0], "q/2" if k % 2 else "(q-2)/2"))

if __name__ == "__main__":
    ks = [int(x) for x in sys.argv[1:]] or [1, 2, 3, 4]
    say("\n=== part4_pure_cubic.py run ===")
    counts()
    for k in ks: run(k)
