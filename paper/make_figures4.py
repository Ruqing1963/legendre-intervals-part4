"""
make_figures4.py -- figures of Part IV (matplotlib, vector PDF output into paper/figures/).

  fig4_phases.pdf : normalised Frobenius eigenvalues zeta_i = exp(2 pi i a_i/120) of the rank-4 Witt sheaf at q = 8,
                    one panel per conjugacy class of pi_lambda (from data/witt_lpoly_k3_b120.txt, all (l3,l5,l1) with
                    multiplicities); the two classes of order 4 are highlighted, the panel titles give S''/(4q^6).
  fig4_kappa.pdf  : the coset parameter kappa_0 of the classes 4 1^2 and 4 2: proportion of pairs with kappa_0 = 0 and 1
                    at q = 8, 32, 128 (from data/part4_kappa_scan.txt).
Usage: python make_figures4.py        (run from paper/; reads ../data/)
"""
import ast
import os
import re
from collections import Counter, defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(os.path.dirname(HERE), "data")
OUT = os.path.join(HERE, "figures")
os.makedirs(OUT, exist_ok=True)

CLASS_NAME = {(1, 1, 1, 1, 1, 1): r"$1^6$", (1, 1, 1, 1, 2): r"$2\,1^4$", (1, 1, 2, 2): r"$2^21^2$", (2, 2, 2): r"$2^3$",
              (1, 1, 1, 3): r"$3\,1^3$", (3, 3): r"$3^2$", (1, 2, 3): r"$3\,2\,1$", (1, 1, 4): r"$4\,1^2$", (2, 4): r"$4\,2$",
              (1, 5): r"$5\,1$", (6,): r"$6$"}
ORDER = [(1, 1, 1, 1, 1, 1), (1, 1, 1, 1, 2), (1, 1, 2, 2), (2, 2, 2), (1, 1, 1, 3), (3, 3), (1, 2, 3), (1, 1, 4), (2, 4),
         (1, 5), (6,)]


def read_phases():
    path = os.path.join(DATA, "witt_lpoly_k3_b120.txt")
    rows = []
    active = False
    for line in open(path, encoding="utf-8"):
        if "normalised Frobenius eigenvalues" in line:
            active = True
            continue
        if active and "small-field Gauss sums" in line:
            break
        m = re.match(r"\s*(\(\(.*\)\))\s+x(\d+)", line)
        if active and m:
            cls, spp, exps = ast.literal_eval(m.group(1))
            rows.append((cls, complex(spp), exps, int(m.group(2))))
    return rows


def fig_phases():
    rows = read_phases()
    by_class = defaultdict(Counter)
    spp_by_class = defaultdict(set)
    for cls, spp, exps, mult in rows:
        for a in exps:
            by_class[cls][a % 120] += mult
        spp_by_class[cls].add(spp.real)
    present = [c for c in ORDER if c in by_class]
    n = len(present)
    cols = 3
    rowsn = int(np.ceil(n / cols))
    fig, axes = plt.subplots(rowsn, cols, figsize=(7.2, 2.4 * rowsn), subplot_kw={"aspect": "equal"})
    axes = np.array(axes).reshape(-1)
    th = np.linspace(0, 2 * np.pi, 400)
    for ax, cls in zip(axes, present):
        ax.plot(np.cos(th), np.sin(th), color="0.75", lw=0.8)
        for a in (0, 30, 60, 90):
            ax.plot([0, np.cos(2 * np.pi * a / 120)], [0, np.sin(2 * np.pi * a / 120)], color="0.9", lw=0.6)
        tot = sum(by_class[cls].values())
        angs = np.array(sorted(by_class[cls]))
        w = np.array([by_class[cls][a] for a in angs], dtype=float) / tot
        hl = cls in ((1, 1, 4), (2, 4))
        ax.scatter(np.cos(2 * np.pi * angs / 120), np.sin(2 * np.pi * angs / 120), s=40 + 600 * w,
                   color="#c0392b" if hl else "#2c3e50", alpha=0.85, edgecolor="white", lw=0.6, zorder=3)
        vals = sorted(spp_by_class[cls])
        lab = ", ".join(f"{v:+g}".replace("+0", "0") if v != 0 else "0" for v in vals)
        ax.set_title(f"{CLASS_NAME[cls]}:  $S''/4q^6\\in\\{{{lab}\\}}$", fontsize=9, color="#c0392b" if hl else "black")
        ax.set_xlim(-1.25, 1.25); ax.set_ylim(-1.25, 1.25)
        ax.set_xticks([]); ax.set_yticks([])
        for s in ax.spines.values():
            s.set_visible(False)
    for ax in axes[n:]:
        ax.set_visible(False)
    fig.suptitle(r"$q=8$: normalised Frobenius eigenvalues $\zeta_i=\exp(2\pi i a_i/120)$ of the Witt sheaf, by class of $\pi_\lambda$"
                 "\n(marker area = multiplicity over all $(\\lambda_3,\\lambda_5,\\lambda_1)$ of the class; order-4 classes in red)",
                 fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(os.path.join(OUT, "fig4_phases.pdf"))
    plt.close(fig)


def read_kappa():
    path = os.path.join(DATA, "part4_kappa_scan.txt")
    res = {}
    q = None
    for line in open(path, encoding="utf-8"):
        m = re.match(r"\s*k=(\d+), q=(\d+)", line)
        if m:
            q = int(m.group(2))
            continue
        m = re.match(r"\s*class (4 1\^2|4 2)\s+(pure|mixed)\s+kappa0=(\d): (\d+)", line)
        if m and q:
            res[(q, m.group(1), m.group(2), int(m.group(3)))] = int(m.group(4))
    return res


def fig_kappa():
    res = read_kappa()
    qs = sorted({k[0] for k in res})
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 2.8), sharey=True)
    for ax, cls, name in zip(axes, ("4 1^2", "4 2"), (r"class $4\,1^2$ (all pairs)", r"class $4\,2$")):
        x = np.arange(len(qs))
        n0 = np.array([res.get((q, cls, "mixed", 0), 0) + res.get((q, cls, "pure", 0), 0) for q in qs], float)
        n1 = np.array([res.get((q, cls, "mixed", 1), 0) + res.get((q, cls, "pure", 1), 0) for q in qs], float)
        tot = n0 + n1
        ax.bar(x - 0.18, n0 / tot, width=0.36, color="#2c3e50", label=r"$\kappa_0=0$")
        ax.bar(x + 0.18, n1 / tot, width=0.36, color="#c0392b", label=r"$\kappa_0=1$")
        for xi, (a, b, t) in enumerate(zip(n0, n1, tot)):
            ax.text(xi - 0.18, a / t + 0.02, f"{int(a)}", ha="center", fontsize=8)
            ax.text(xi + 0.18, b / t + 0.02, f"{int(b)}", ha="center", fontsize=8)
        ax.set_xticks(x); ax.set_xticklabels([f"$q={q}$" for q in qs])
        ax.set_title(name, fontsize=10)
        ax.set_ylim(0, 1.15)
        ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    axes[0].set_ylabel("proportion of pairs")
    axes[1].legend(frameon=False, fontsize=9, loc="upper right")
    fig.suptitle(r"The coset parameter $\kappa_0=\frac{1}{2}\mathfrak{q}^{(4)}_{\lambda_1=0}(r_1)$: $s_4(\lambda_1)\neq0$ iff $\mathrm{Tr}(\lambda_1\tau_0)=\kappa_0$"
                 "\n(numbers = pairs $(\\lambda_3,\\lambda_5)$; $\\kappa_0$ is not constant on the class for $q \\geq 32$)", fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 0.88))
    fig.savefig(os.path.join(OUT, "fig4_kappa.pdf"))
    plt.close(fig)


if __name__ == "__main__":
    fig_phases()
    fig_kappa()
    print("written:", sorted(os.listdir(OUT)))
