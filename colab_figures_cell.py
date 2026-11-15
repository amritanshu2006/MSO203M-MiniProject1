# @title Extra step: make the report figures (run after Steps 1-5)
# Uses only the tools loaded in Step 1 (make_g, target, exact_solution, evolve)
# and the case from Step 2. Saves the 5 PNG files used in the report.
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from PIL import Image

c = case                      # assigned data from Step 2 (roll 250120)
T = c.T; eps = c.epsilon
out = Path("report_figures"); out.mkdir(exist_ok=True)

# Boundary functions (same formulas as Steps 3 and 4)
F1 = "(2*exp(-t))^kappa*(pos(0.65*bump(2*exp(-t)-1,0.25,0.15)-0.26)+0.72*0.9*bump(2*exp(-t)-1,xstar,w))"
F2 = "(2*exp(-t))^kappa*(0.6*0.65*bump(2*exp(-t)-1,0.25,0.15)+0.72*0.9*bump(2*exp(-t)-1,xstar,w))"
FS = "(2*exp(-t))^kappa*pos(D(2*exp(-t)-1)-0.25)"   # alternative (ii): shift and clip
FA = "0.61*(2*exp(-t))^kappa*D(2*exp(-t)-1)"        # rejected (i): uniform scaling
FD = "(2*exp(-t))^kappa*D(2*exp(-t)-1)"             # g_D: exact reproduction
g1 = make_g(F1, c, "original"); g2 = make_g(F2, c, "narrow")
gS = make_g(FS, c, "original"); gA = make_g(FA, c, "original"); gD = make_g(FD, c, "original")

PURPLE = "#6a3d9a"; ORANGE = "#e66101"; TEAL = "#1b9e77"; GREY = "#777777"
LIGHT = "#ece4f3"   # light purple fill

# ---------- Figure 1: characteristics sketch ----------
plt.rcParams.update(plt.rcParamsDefault)
t0 = 0.45
fig, ax = plt.subplots(figsize=(6.4, 4.6), dpi=200)
tt = np.linspace(0, T, 400); sep = np.exp(tt) - 1
ax.fill_betweenx(tt, 0, sep, color=LIGHT)
ax.fill_betweenx(tt, sep, 1, color="#f3f3f3")
for s in np.arange(0.1, T, 0.1):                     # characteristics from the boundary (0,s)
    t = np.linspace(s, T, 200); ax.plot(np.exp(t - s) - 1, t, color=PURPLE, lw=1)
for s in np.arange(0.1, 1.0, 0.1):                   # characteristics from the initial line (s,0)
    t = np.linspace(0, T, 200); xs = (1 + s) * np.exp(t) - 1; m = xs <= 1
    ax.plot(xs[m], t[m], color="#888888", lw=1)
ax.plot(sep, tt, color="k", lw=2.2)                  # separating curve x = e^t - 1
ax.axhline(t0, color=ORANGE, lw=1.2, ls="--")
xc = np.exp(t0) - 1
ax.plot([xc], [t0], "o", color=ORANGE, ms=5)
ax.annotate(r"$t_0=0.45$: $x_c=e^{0.45}-1\approx0.568$", (xc, t0), (0.60, 0.36), fontsize=8.5,
            bbox=dict(fc="white", ec="none", alpha=0.9, pad=1),
            arrowprops=dict(arrowstyle="-", lw=0.7, color=ORANGE))
ax.axhline(T, color="k", lw=0.8)
ax.text(0.01, T + 0.008, r"$T=\log 2$", fontsize=9)
ax.text(0.05, 0.55, "boundary region\n$x<e^t-1$\n$u=(1+x)^{-\\kappa}g(t-\\log(1+x))$", fontsize=8.5,
        color="#4a2a6e", bbox=dict(fc="white", ec=PURPLE, lw=0.6, alpha=0.95))
ax.text(0.55, 0.08, "initial region $x\\geq e^t-1$\n$u=0$", fontsize=8.5, color="#444",
        bbox=dict(fc="white", ec="#888", lw=0.6, alpha=0.95))
ax.text(0.30, 0.20, r"$x=e^t-1$ (from $(0,0)$)", fontsize=8.5, rotation=38,
        bbox=dict(fc="white", ec="none", alpha=0.9, pad=1))
ax.set_xlim(0, 1); ax.set_ylim(0, T + 0.05)
ax.set_xlabel("x"); ax.set_ylabel("t")
ax.set_title(r"Characteristics of $u_t+(1+x)u_x=-\kappa u$", fontsize=10)
for sp in ("top", "right"): ax.spines[sp].set_visible(False)
fig.tight_layout(); fig.savefig(out / "partA_characteristics.png"); plt.show()

# ---------- Figure 2: boundary data and tolerance band ----------
plt.rcParams.update({"font.size": 8.5, "axes.spines.top": False, "axes.spines.right": False})
x = np.linspace(0, 1, 40001); tt = np.linspace(0, T, 40001)
fig, ax = plt.subplots(1, 2, figsize=(7.4, 2.55), dpi=220)
a = ax[0]
a.fill_between(tt, 0, gD(tt), color="#e9e9e9", label=r"$g_D$ (exact reproduction), area $\approx0.152$")
a.plot(tt, gD(tt), color=GREY, lw=1)
a.plot(tt, g1(tt), color=PURPLE, lw=1.8, label=r"chosen $g$, area $\approx0.0837$")
a.plot(tt, g2(tt), color=TEAL, lw=1.1, ls="--", label=r"$g_n$ (narrow target), area $\approx0.0649$")
a.axhline(1, color=ORANGE, ls="--", lw=1); a.text(0.26, 0.96, r"cap $g\leq1$", color=ORANGE, va="top")
i = np.argmax(gD(tt))
a.annotate(r"$g_D(s^*)=1.229>1$", (tt[i], gD(tt)[i]), (0.14, 1.17), arrowprops=dict(arrowstyle="-", lw=.6))
a.set_xlabel("boundary time $t$"); a.set_ylabel("$g(t)$"); a.set_ylim(0, 2.05); a.set_xlim(0, T)
a.legend(loc="upper right", fontsize=7, frameon=False); a.set_title("(a) boundary data", fontsize=9)
b = ax[1]
D = target(c, x, "original")
b.fill_between(x, np.maximum(D - eps, 0), D + eps, color=LIGHT, label=r"allowed band $D\pm\varepsilon$")
b.plot(x, D, color="k", lw=1.2, label="target $D$")
b.plot(x, exact_solution(c, gA, x, T), color=ORANGE, lw=1, ls=":", label=r"rejected: $0.61\,D$ ($E=0.351$)")
b.plot(x, exact_solution(c, g1, x, T), color=PURPLE, lw=1.8, label=r"chosen: $u_g(x,T)=H(x)$ ($E=0.26$)")
b.set_xlabel("$x$"); b.set_ylabel("value at $t=T$"); b.set_xlim(0, 1); b.set_ylim(0, 1.25)
b.legend(loc="upper left", fontsize=6.8, frameon=False)
b.set_title("(b) exact final profile vs. tolerance band", fontsize=9)
fig.tight_layout(); fig.savefig(out / "figB_design.png"); plt.show()

# ---------- Figure 4: delta_N vs N and zoom on the second peak ----------
Ns = [50, 100, 200, 400, 800, 1600]
res = {}
for name, g, which in (("original, chosen g", g1, "original"),
                       ("original, shift-clip alternative", gS, "original"),
                       ("narrow, chosen $g_n$", g2, "narrow")):
    ue = exact_solution(c, g, x, T); d = []
    for N in Ns:
        xx, U, M, dt = evolve(c, g, N)
        d.append(float(np.abs(np.interp(x, xx, U) - ue).max()))
    res[name] = d
    print(name, [round(v, 4) for v in d])
fig, ax = plt.subplots(1, 2, figsize=(7.4, 2.5), dpi=220)
a = ax[0]
for (name, d), col, mk in zip(res.items(), (PURPLE, TEAL, ORANGE), ("o", "s", "^")):
    a.loglog(Ns, d, marker=mk, color=col, lw=1.2, ms=4, label=name)
NN = np.array(Ns, float)
a.loglog(NN, 0.2 * (NN / 400) ** -1, color="k", lw=.7, ls="--", label="slope $-1$ (first order)")
a.loglog(NN, 0.2 * (NN / 400) ** -0.5, color="k", lw=.7, ls=":", label="slope $-1/2$")
for N in (100, 200, 400): a.axvline(N, color="#cccccc", lw=.6, zorder=0)
a.set_xlabel("$N$ (spatial intervals)"); a.set_ylabel(r"$\delta_N$")
a.legend(fontsize=6.3, frameon=False, loc="lower left")
a.set_title(r"(a) numerical error $\delta_N$ vs. $N$", fontsize=9)
b = ax[1]
xz = np.linspace(0.68, 0.95, 3000)
for which, g, col, lab in (("original", g1, PURPLE, "original"), ("narrow", g2, ORANGE, "narrow")):
    b.plot(xz, exact_solution(c, g, xz, T), color=col, lw=1.6, label=f"exact, {lab}")
    xx, U, M, dt = evolve(c, g, 400); m = (xx >= 0.68) & (xx <= 0.95)
    b.plot(xx[m], U[m], color=col, lw=.9, ls="--", marker=".", ms=2.5, label=f"$N=400$, {lab}")
b.set_xlabel("$x$"); b.set_ylabel("$u(x,T)$")
b.set_title(r"(b) zoom on the second peak ($N=400$ nodes)", fontsize=9)
b.legend(fontsize=6.5, frameon=False, loc="upper right")
fig.tight_layout(); fig.savefig(out / "figC_refinement.png"); plt.show()
plt.rcParams.update(plt.rcParamsDefault)

# ---------- Figure 3: right halves of the notebook's own Step 3/4 plots ----------
for n in ("original", "narrow"):
    im = Image.open(Path(results_folder) / f"{n}_comparison.png")
    w, h = im.size
    im.crop((int(w * 0.5), 0, w, h)).save(out / f"{n}_right.png")

print("Saved:", sorted(p.name for p in out.glob("*.png")))
try:
    from google.colab import files
    import shutil
    shutil.make_archive("report_figures", "zip", out)
    files.download("report_figures.zip")
except ImportError:
    print("Files are in the folder", out)
