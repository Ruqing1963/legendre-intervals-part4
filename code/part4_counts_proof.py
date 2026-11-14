r"""Verification of the proof of the class counts (Target Theorem 1.1) from Lemma 3.14, the first moment (Part III,
Theorem 1.3) and the third moment (Theorem 2.3).

Incidences (sextic, Galois-stable pair of roots {r5, r6}) <-> (y, s) in (F_q \ {0,1}) x F_q^x (Lemma 3.14):
  tau0 = y^2 + y, rho0 = s tau0^2, z = y^4, b = rho0 z / tau0, a = b^2/rho0 + tau0^4 + rho0 tau0^2 + rho0^2,
  cofactor K(u) = (u^6 + a u^2 + b u + b^2)/(u^2 + tau0 u + rho0) = u^4 + tau0 u^3 + (tau0^2 + rho0) u^2 + tau0^3 u + b^2/rho0,
  depressed quartic Phi(x) = x^4 + a' x^2 + b' x + c', a' = 1/(1+t), b' = c' = a'/s, t = (y/(y+1))^4,
  resolvent Y^3 + a' Y + b'.  A rational root r gives u = r s, xi = 1/r, with
  Tr(s) = Tr(u xi) (splitting of the quadratic), Tr(u/(u+1)) = 0 iff r is the only rational root,
  Tr((u+1) xi) = Tr(c'/(a'+r^2)^2) = 0 iff the partition quadratic X^2 + (a'+r^2) X + c' splits.
The script tallies the incidences (y, s, r) by (number of resolvent roots, Tr(s), Tr((u+1)xi), class of the sextic),
compares with the predicted cell counts, classifies every sextic directly (gcd with x^{q^m} - x, m = 1,2,3), and
checks the polynomial identities with sympy.   Usage: python part4_counts_proof.py [k ...]   (default 3 5 7)."""
import sys, os, time, itertools
import numpy as np
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = open(os.path.join(HERE, "..", "data", "part4_counts_proof.txt"), "a", encoding="utf8")
def say(s=""):
    print(s, flush=True); OUT.write(s + "\n"); OUT.flush()

# ---------------- GF(2^k) ----------------
IRRED = {3: 0b1011, 5: 0b100101, 7: 0b10000011, 9: 0b1000010001}
class GF:
    def __init__(self, k):
        self.k = k; self.q = 1 << k; self.mod = IRRED[k]
        q = self.q
        def clmul(a, b):
            r = 0
            while b:
                if b & 1: r ^= a
                a <<= 1
                if a & q: a ^= self.mod
                b >>= 1
            return r
        # primitive element
        for g in range(2, q):
            x, n = 1, 0
            seen = set()
            while True:
                x = clmul(x, g); n += 1
                if x == 1: break
            if n == q - 1: break
        self.exp = np.zeros(2 * q, np.int64); self.log = np.zeros(q, np.int64)
        x = 1
        for i in range(q - 1):
            self.exp[i] = x; self.log[x] = i; x = clmul(x, g)
        self.exp[q - 1:2 * q - 2] = self.exp[:q - 1]
        self.tr = np.zeros(q, np.int64)
        for x in range(q):
            t, y = 0, x
            for _ in range(k):
                t ^= y; y = clmul(y, y)
            self.tr[x] = t
        assert set(self.tr) == {0, 1}
    def mul(self, a, b):
        if a == 0 or b == 0: return 0
        return int(self.exp[self.log[a] + self.log[b]])
    def inv(self, a): return int(self.exp[(self.q - 1 - self.log[a]) % (self.q - 1)])
    def div(self, a, b): return self.mul(a, self.inv(b))
    def pw(self, a, n):
        if a == 0: return 0 if n else 1
        return int(self.exp[(self.log[a] * n) % (self.q - 1)])
    def sqrt(self, a): return self.pw(a, self.q // 2)
    # vectorised: product of arrays (elements), mul by scalar
    def vmul(self, A, B):
        A = np.asarray(A); B = np.asarray(B)
        r = self.exp[self.log[A] + self.log[B]]
        return np.where((A == 0) | (B == 0), 0, r)

# ---------------- polynomials over GF (low degree first, lists of ints) ----------------
def pmul(F, a, b):
    r = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                if y: r[i + j] ^= F.mul(x, y)
    return r
def pmod(F, a, m):
    a = a[:]
    dm = len(m) - 1; lead = F.inv(m[-1])
    while len(a) - 1 >= dm:
        c = a[-1]
        if c:
            f = F.mul(c, lead)
            for i in range(dm + 1):
                a[len(a) - 1 - dm + i] ^= F.mul(f, m[i])
        a.pop()
        while len(a) > 1 and a[-1] == 0: a.pop()
        if len(a) == 1 and a[0] == 0: return [0]
        if len(a) == 0: return [0]
    return a if a else [0]
def pdivmod(F, a, m):
    a = a[:]; dm = len(m) - 1; lead = F.inv(m[-1]); qt = [0] * max(1, len(a) - dm)
    for i in range(len(a) - 1, dm - 1, -1):
        c = a[i]
        if c:
            f = F.mul(c, lead); qt[i - dm] = f
            for j in range(dm + 1): a[i - dm + j] ^= F.mul(f, m[j])
    r = a[:dm]
    while len(r) > 1 and r[-1] == 0: r.pop()
    return qt, r
def pgcd(F, a, b):
    while not (len(b) == 1 and b[0] == 0):
        a, b = b, pmod(F, a, b)
    return a
def xq_pow_mod(F, m, e):
    """x^(2^e) mod m by repeated squaring"""
    x = [0, 1]
    if len(m) <= 2: return pmod(F, x, m)
    for _ in range(e): x = pmod(F, pmul(F, x, x), m)
    return x
def gcd_degrees(F, f, ms):
    """degrees of gcd(f, x^{q^m} - x) for m in ms"""
    out = []
    for m in ms:
        xp = xq_pow_mod(F, f, F.k * m)
        d = xp[:]
        if len(d) < 2: d = d + [0] * (2 - len(d))
        d[1] ^= 1
        while len(d) > 1 and d[-1] == 0: d.pop()
        g = pgcd(F, f, d)
        out.append(len(g) - 1)
    return out

SIG = {(6, 6, 6): "1^6", (4, 6, 4): "2 1^4", (2, 6, 2): "2^2 1^2", (0, 6, 0): "2^3", (3, 3, 6): "3 1^3", (0, 0, 6): "3^2",
       (1, 3, 4): "3 2 1", (2, 2, 2): "4 1^2", (0, 2, 0): "4 2", (1, 1, 1): "5 1", (0, 0, 0): "6"}
QSIG = {(4, 4): "1111", (2, 4): "211", (0, 4): "22", (1, 1): "13", (0, 0): "4"}
CLASSES = ["1^6", "2 1^4", "2^2 1^2", "2^3", "3 1^3", "3^2", "3 2 1", "4 1^2", "4 2", "5 1", "6"]
STAB = {"1^6": 15, "2 1^4": 7, "2^2 1^2": 3, "2^3": 3, "3 1^3": 3, "3^2": 0, "3 2 1": 1, "4 1^2": 1, "4 2": 1, "5 1": 0, "6": 0}
Q = sp.symbols("q")
TABLE = {"1^6": (Q - 2) * (Q - 8) / 720, "2 1^4": (Q - 2) * (Q - 8) / 48, "2^2 1^2": Q * (Q - 2) / 16, "2^3": (Q - 2) * (Q - 8) / 48 + (Q - 2) / 6,
         "3 1^3": (Q + 1) * (Q - 2) / 18, "3^2": (Q + 1) * (Q - 2) / 18, "3 2 1": (Q + 1) * (Q - 2) / 6, "4 1^2": (Q - 2) * (Q - 4) / 8 + (Q - 2) / 2,
         "4 2": Q * (Q - 2) / 8, "5 1": (Q ** 2 + 1) / 5, "6": (Q + 1) * (Q - 2) / 6 + (Q + 1) / 3}       # totals (mixed + pure)
F1 = {"1^6": 6, "2 1^4": 4, "2^2 1^2": 2, "2^3": 0, "3 1^3": 3, "3^2": 0, "3 2 1": 1, "4 1^2": 2, "4 2": 0, "5 1": 1, "6": 0}
F3 = {"1^6": 6, "2 1^4": 4, "2^2 1^2": 2, "2^3": 0, "3 1^3": 6, "3^2": 6, "3 2 1": 4, "4 1^2": 2, "4 2": 0, "5 1": 1, "6": 0}

# predicted cell counts: key (number of resolvent roots, Tr(s), Tr((u+1)xi)) -> {class: multiplicity}
CELLS = {(1, 0, 1): {"4 1^2": 1}, (1, 1, 1): {"4 2": 1}, (1, 0, 0): {"2 1^4": 6}, (1, 1, 0): {"2^2 1^2": 2},
         (3, 1, 1): {"2^3": 6}, (3, 0, 0): {"1^6": 45, "2^2 1^2": 1}, (3, 1, 0): {"2 1^4": 3, "2^3": 3}, (3, 0, 1): {"2^2 1^2": 2}}
CELLSIZE = {(1, 0, 1): (Q - 2) / 2 * Q / 4, (1, 1, 1): (Q - 2) / 2 * Q / 4, (1, 0, 0): (Q - 2) / 2 * (Q / 4 - 2), (1, 1, 0): (Q - 2) / 2 * Q / 4,
            (3, 1, 1): (Q - 2) / 2 * Q / 4, (3, 0, 0): (Q - 2) / 2 * (Q / 4 - 1), (3, 1, 0): (Q - 2) / 2 * (Q / 4 - 1), (3, 0, 1): (Q - 2) / 2 * Q / 4}

def run(k):
    t0 = time.time(); F = GF(k); q = F.q
    say("\n=== k = %d, q = %d ===" % (k, q))
    # ---- direct classification of all sextics (skipped for q > 128, where Part III did it; the table values are used) ----
    counts = {c: 0 for c in CLASSES}; cls_of = {}
    if q <= 128:
        for a in range(q):
            for b in range(1, q):
                f = [F.mul(b, b), b, a, 0, 0, 0, 1]
                sig = tuple(gcd_degrees(F, f, (1, 2, 3)))
                c = SIG[sig]; counts[c] += 1; cls_of[(a, b)] = c
        ok = all(counts[c] == int(TABLE[c].subs(Q, q)) for c in CLASSES)
        say("direct classification of all %d sextics: " % (q * (q - 1)) + ", ".join("%s:%d" % (c, counts[c]) for c in CLASSES))
        say("  equals Table (totals, mixed + pure) for every class: %s   [%.0f s]" % (ok, time.time() - t0))
    else:
        counts = {c: int(TABLE[c].subs(Q, q)) for c in CLASSES}
        say("direct classification skipped (Part III, Proposition 5.2, q <= 512); table values used for the multiplicities")
    inc_pred = sum(counts[c] * STAB[c] for c in CLASSES)
    say("  incidences (sextic, stable pair) from the classes: %d; (q-1)(q-2) = %d: %s" % (inc_pred, (q - 1) * (q - 2), inc_pred == (q - 1) * (q - 2)))
    # ---- the (y, s) parametrisation ----
    seen = set(); cells = {}; nroots_hist = {}; bad = 0; typecheck = 0
    quad_u = {}
    one = 1
    for y in range(2, q):
        tau0 = F.mul(y, y) ^ y; z = F.pw(y, 4)
        t = F.pw(F.div(y, y ^ 1), 4); ap = F.inv(t ^ 1)
        for s in range(1, q):
            rho0 = F.mul(s, F.mul(tau0, tau0)); b = F.div(F.mul(rho0, z), tau0)
            b2 = F.mul(b, b)
            a = F.div(b2, rho0) ^ F.pw(tau0, 4) ^ F.mul(rho0, F.mul(tau0, tau0)) ^ F.mul(rho0, rho0)
            f = [b2, b, a, 0, 0, 0, 1]
            qd = [rho0, tau0, 1]
            K, rem = pdivmod(F, f, qd)
            assert rem == [0], "quadratic does not divide"
            Kl = [F.div(b2, rho0), F.pw(tau0, 3), F.mul(tau0, tau0) ^ rho0, tau0, 1]
            assert K == Kl, "cofactor formula"
            key = (a, b, tau0, rho0); assert key not in seen, "incidence repeated"; seen.add(key)
            # Phi and its type; the type of K must agree
            bp = F.div(ap, s); Phi = [bp, bp, ap, 0, 1]
            tP = QSIG[tuple(gcd_degrees(F, Phi, (1, 2)))]; tK = QSIG[tuple(gcd_degrees(F, K, (1, 2)))]
            typecheck += (tP == tK)
            # class from the quadratic splitting (Tr(s) = 0 iff splits) and the type of Phi
            spl = F.tr[s] == 0
            pred = {"1111": "1^6" if spl else "2 1^4", "211": "2 1^4" if spl else "2^2 1^2", "22": "2^2 1^2" if spl else "2^3",
                    "13": "3 1^3" if spl else "3 2 1", "4": "4 1^2" if spl else "4 2"}[tP]
            c = cls_of.get((a, b), pred)
            assert pred == c, "class mismatch"
            # rational roots of the resolvent Y^3 + a' Y + b'
            roots = [r for r in range(q) if F.mul(r, F.mul(r, r)) ^ F.mul(ap, r) ^ bp == 0]
            nr = len(roots); nroots_hist[(tP, nr)] = nroots_hist.get((tP, nr), 0) + 1
            for r in roots:
                assert r != 0
                u = F.mul(r, s); xi = F.inv(r)
                assert u not in (0, 1)
                xi0 = F.sqrt(F.div(u, u ^ 1)); assert xi != xi0
                assert F.tr[F.mul(u, xi)] == F.tr[s]
                assert F.tr[F.div(u, u ^ 1)] == (0 if nr == 1 else 1)
                A2 = ap ^ F.mul(r, r)
                assert F.tr[F.div(bp, F.mul(A2, A2))] == F.tr[F.mul(u ^ 1, xi)]
                key2 = (nr, int(F.tr[s]), int(F.tr[F.mul(u ^ 1, xi)]))
                cells.setdefault(key2, {}); cells[key2][c] = cells[key2].get(c, 0) + 1
                quad_u.setdefault(u, set()).add((xi, key2[1], key2[2]))
    say("(y, s) plane: %d incidences, all distinct, quadratic divides and cofactor agrees with Lemma 3.14 in every case; type(Phi) = type(K) in %d cases"
        % (len(seen), typecheck))
    say("  resolvent roots by type of Phi: " + ", ".join("%s:%d roots x%d" % (tp, nr, n) for (tp, nr), n in sorted(nroots_hist.items())))
    allok = True
    for key2 in sorted(CELLS):
        got = cells.get(key2, {}); pred = CELLS[key2]; size = int(CELLSIZE[key2].subs(Q, q))
        tot = sum(got.values())
        exp = {c: pred[c] * counts[c] for c in pred if pred[c] * counts[c]}
        okc = (got == exp) and (tot == size)
        allok &= okc
        say("  cell (roots=%d, Tr s=%d, Tr((u+1)xi)=%d): %s; size %d (predicted %d) %s" % (key2 + (dict(sorted(got.items())), tot, size, "OK" if okc else "MISMATCH")))
    say("  all eight cells as predicted: %s" % allok)
    # quadrant distribution per u
    qok = True
    for u, S in quad_u.items():
        cnt = {}
        for (xi, t1, t2) in S: cnt[(t1, t2)] = cnt.get((t1, t2), 0) + 1
        tr_u = F.tr[F.div(u, u ^ 1)]
        exp = {(0, 0): q // 4 - 2, (0, 1): q // 4, (1, 0): q // 4, (1, 1): q // 4} if tr_u == 0 else {(0, 0): q // 4 - 1, (0, 1): q // 4, (1, 0): q // 4 - 1, (1, 1): q // 4}
        qok &= (cnt == {kk: v for kk, v in exp.items() if v})
    say("  (u, xi) bijection: %d values of u (= q-2: %s), quadrant counts per u as predicted: %s" % (len(quad_u), len(quad_u) == q - 2, qok))
    # the deductions
    m1 = sum(F1[c] * counts[c] for c in CLASSES); m3 = sum(F3[c] * counts[c] for c in CLASSES)
    say("  first moment %d = (q-1)^2 = %d: %s;  third moment %d = 2q^2-3q-1 = %d: %s" % (m1, (q - 1) ** 2, m1 == (q - 1) ** 2, m3, 2 * q * q - 3 * q - 1, m3 == 2 * q * q - 3 * q - 1))
    say("  [%.0f s]" % (time.time() - t0))

def symbolic():
    say("=== polynomial identities (sympy) ===")
    N = {c: TABLE[c] for c in CLASSES}
    say("  sum N = q(q-1): %s" % (sp.simplify(sum(N.values()) - Q * (Q - 1)) == 0))
    say("  sum N f1 = (q-1)^2: %s" % (sp.simplify(sum(F1[c] * N[c] for c in CLASSES) - (Q - 1) ** 2) == 0))
    say("  sum N f3 = 2q^2-3q-1: %s" % (sp.simplify(sum(F3[c] * N[c] for c in CLASSES) - (2 * Q ** 2 - 3 * Q - 1)) == 0))
    say("  sum N stab = (q-1)(q-2): %s" % (sp.simplify(sum(STAB[c] * N[c] for c in CLASSES) - (Q - 1) * (Q - 2)) == 0))
    for key2, pred in CELLS.items():
        say("  cell %s: sum mult*N = size: %s" % (str(key2), sp.simplify(sum(m * N[c] for c, m in pred.items()) - CELLSIZE[key2]) == 0))
    # the two slices
    s0 = 15 * N["1^6"] + 6 * N["2 1^4"] + N["2^2 1^2"] + N["4 1^2"] + 3 * N["3 1^3"]
    s1 = N["2 1^4"] + 2 * N["2^2 1^2"] + 3 * N["2^3"] + N["4 2"] + N["3 2 1"]
    say("  slice Tr s = 0: (q-2)^2/2: %s;  slice Tr s = 1: q(q-2)/2: %s" % (sp.simplify(s0 - (Q - 2) ** 2 / 2) == 0, sp.simplify(s1 - Q * (Q - 2) / 2) == 0))
    # solving in the order of the proof
    n412 = (Q - 2) / 2 * Q / 4; n42 = n412; n214 = (Q - 2) / 2 * (Q / 4 - 2) / 6; n2212 = (Q - 2) / 2 * Q / 4 / 2
    n23 = (Q - 2) / 2 * Q / 4 / 6; n16 = ((Q - 2) / 2 * (Q / 4 - 1) - n2212) / 45
    n313 = ((Q - 2) ** 2 / 2 - (15 * n16 + 6 * n214 + n2212 + n412)) / 3
    n321 = Q * (Q - 2) / 2 - (n214 + 2 * n2212 + 3 * n23 + n42)
    n51 = (Q - 1) ** 2 - (6 * n16 + 4 * n214 + 2 * n2212 + 3 * n313 + n321 + 2 * n412)
    n32 = (2 * Q ** 2 - 3 * Q - 1 - (6 * n16 + 4 * n214 + 2 * n2212 + 6 * n313 + 4 * n321 + 2 * n412 + n51)) / 6
    n6 = Q * (Q - 1) - (n16 + n214 + n2212 + n23 + n313 + n32 + n321 + n412 + n42 + n51)
    sol = {"1^6": n16, "2 1^4": n214, "2^2 1^2": n2212, "2^3": n23, "3 1^3": n313, "3^2": n32, "3 2 1": n321, "4 1^2": n412, "4 2": n42, "5 1": n51, "6": n6}
    for c in CLASSES:
        say("  solved N(%s) = %s  equals table: %s" % (c, sp.factor(sol[c]), sp.simplify(sol[c] - N[c]) == 0))

if __name__ == "__main__":
    ks = [int(x) for x in sys.argv[1:]] or [3, 5, 7]
    say("\n=== part4_counts_proof.py run ===")
    symbolic()
    for k in ks: run(k)
