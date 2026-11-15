# Boundary Data for a First-Order PDE

**MSO203M: Introduction to Partial Differential Equations, IIT Kanpur | Mini-project 1**
**Author:** Amritanshu Kumar (Roll No. 250120) | **Case code:** 16931

In this project I choose the boundary input of a first-order transport PDE so that the solution at a final time matches a target profile within a given tolerance. I solve the PDE exactly with the method of characteristics, prove that exact matching is impossible for my assigned case, design an admissible boundary function with a proven error bound, and then study how well the supplied upwind finite-difference scheme reproduces the exact solution.

---

## Table of contents
1. [Problem statement](#1-problem-statement)
2. [My assigned data](#2-my-assigned-data)
3. [Part A: Characteristic solution](#3-part-a-characteristic-solution)
4. [Part B: Choosing the boundary data](#4-part-b-choosing-the-boundary-data)
5. [Part C: The numerical method](#5-part-c-the-numerical-method)
6. [Part D: Narrower target](#6-part-d-narrower-target)
7. [Conclusions](#7-conclusions)
8. [Figures](#8-figures)
9. [How to reproduce the results](#9-how-to-reproduce-the-results)
10. [Repository structure](#10-repository-structure)
11. [Tools used](#11-tools-used)
12. [Acknowledgements](#12-acknowledgements)

---

## 1. Problem statement

For a boundary function $g$, let $u_g$ solve

$$
u_t + (1+x)\,u_x = -\kappa u,\qquad 0<x<1,\; 0<t\le T=\log 2,
$$
$$
u(x,0)=0,\qquad u(0,t)=g(t).
$$

The target at the final time is a two-peak profile

$$
D(x) = 0.65\,b_{0.25,\,0.15}(x) + 0.90\,b_{x^*,\,w}(x),\qquad
b_{a,h}(x)=\begin{cases}\cos^2\!\dfrac{\pi(x-a)}{2h}, & |x-a|<h,\\ 0, & \text{otherwise.}\end{cases}
$$

The boundary function must be **admissible**, meaning:

- $g$ is continuous and piecewise $C^1$ (corners are allowed, jumps are not);
- $g(0)=g(T)=0$ and $0\le g\le 1$;
- $\int_0^T g(t)\,dt \le B$ (the budget).

It must also meet the accuracy requirement

$$
E(g)=\max_{0\le x\le 1}\,|u_g(x,T)-D(x)| \le \varepsilon .
$$

The project has four parts:

| Part | Task |
|---|---|
| A | Find the characteristic solution for an arbitrary $g$ |
| B | Show whether exact reproduction is possible, then construct an acceptable $g$ |
| C | Analyse the supplied numerical method and run it for $N=100,200,400$ |
| D | Repeat the design and the computation for a narrower second peak |

---

## 2. My assigned data

The course notebook generates the parameters from the roll number:

| Symbol | Meaning | Value |
|---|---|---|
| $\kappa$ | damping coefficient | 0.525 |
| $B$ | budget on $\int g$ | 0.093 |
| $x^*$ | centre of the second peak | 0.81 |
| $w$ | half-width of the second peak | 0.095 |
| $w_n$ | narrower half-width (Part D) | 0.025 |
| $\varepsilon$ | tolerance on $E(g)$ | 0.29 |
| $T$ | final time | $\log 2 \approx 0.693147$ |
| $R^*=(1+x^*)^{-\kappa}$ | largest value reachable at $x^*$ | 0.7323501 |

---

## 3. Part A: Characteristic solution

### A1. Characteristic equations and solution formula

Along a characteristic $(X(r,s),T(r,s))$ with $V=u(X,T)$, where $s$ labels the starting point and $r$ is the elapsed time:

$$
X_r = 1+X,\qquad T_r = 1,\qquad V_r=-\kappa V
\;\;\Longrightarrow\;\;
1+X=(1+X_0)e^{r},\quad T=T_0+r,\quad V=V_0e^{-\kappa r}.
$$

There are two families:

- **From the initial line** $(s,0)$: $X=(1+s)e^r-1$, $T=r$, $V=0$. This family covers the points with $x \ge e^t-1$.
- **From the boundary** $(0,s)$: $X=e^r-1$, $T=s+r$, $V=e^{-\kappa r}g(s)$, with travel time $r=\log(1+x)$. This family covers the points with $x < e^t-1$.

Together they give

$$
\boxed{\,u_g(x,t)=\begin{cases}0, & x\ge e^t-1 \quad\text{(initial data)},\\[2pt]
(1+x)^{-\kappa}\,g\big(t-\log(1+x)\big), & x<e^t-1 \quad\text{(boundary data)}.\end{cases}}
$$

- The two regions are separated by the corner characteristic $x=e^t-1$. Both formulas give 0 on it because $g(0)=0$, so $u_g$ is continuous.
- **No condition is needed at $x=1$.** The speed $1+x>0$, so every characteristic leaves through $x=1$ and none enters there. Each point traces back to $t=0$ or to $x=0$.

### A2. Dependence on $g$

- At $t_0=0.45$, only $0\le x<e^{0.45}-1\approx 0.568$ can depend on $g$. For larger $x$, $u=0$.
- At $T=\log 2$, $e^T-1=1$, so **every interior point depends on $g$**. The final profile is

$$
u_g(x,T)=(1+x)^{-\kappa}\,g\!\left(\log\frac{2}{1+x}\right).
$$

Each final position reads $g$ at a single time. The peak at $x^*=0.81$ is set by $g$ near $s^*=\log(2/1.81)=0.0998$, and the peak at $x=0.25$ by $g$ near $s=0.470$.

---

## 4. Part B: Choosing the boundary data

### B1. Exact reproduction is impossible

Inverting the final-profile formula gives the only $g$ that reproduces $D$ exactly:

$$
g_D(s)=(2e^{-s})^{\kappa}\,D(2e^{-s}-1).
$$

This $g_D$ fails two restrictions:

- **The cap $g\le1$.** $g_D(s^*)=0.9/R^*=1.229>1$. In fact, for **every** admissible $g$, $u_g(x^*,T)\le R^*=0.7324<0.9$, so $E(g)\ge 0.1676$. The damping along the characteristic means even $g=1$ cannot reach the peak.
- **The budget.** Changing variables gives $\int_0^T g\,ds=\int_0^1 (1+x)^{\kappa-1}H(x)\,dx$, where $H=u_g(\cdot,T)$. A lower bound then gives $\int g_D \ge 0.132 > B$ (the numerical value is $0.1523$).

### B2. Construction of $g$

**Approach.** I design the final profile $H(x)=u_g(x,T)$ first and recover $g$ from

$$
g(s)=(2e^{-s})^{\kappa}H(2e^{-s}-1).
$$

All the restrictions on $g$ then become conditions on $H$.

I compared three candidate profiles:

| Candidate $H$ | $E(g)$ | $\max g$ | $\int g$ | Result |
|---|---|---|---|---|
| (i) uniform scaling $\alpha D$ | $0.9(1-\alpha)$ | $1.229\alpha$ | $0.152\alpha$ | **Fails**: the tolerance needs $\alpha\ge0.678$, but the budget needs $\alpha\le0.611$ |
| (ii) shift and clip $(D-0.25)_+$ | 0.25 | 0.888 | 0.0757 | Valid, but narrows both peaks ($\delta_{400}=0.255$) |
| (iii) **chosen**: clip peak 1, scale peak 2 | **0.26** | **0.885** | **0.0837** | Valid, and peak 2 keeps its full width ($\delta_{400}=0.199$) |

**Why option (iii).** The cap affects only peak 2, and peak 2 is also the hardest for the numerical scheme. So I keep peak 2 smooth and full-width by scaling it by 0.72. To stay within the budget, I trim mainly peak 1, where each unit of height costs the most area. With $\beta_1=0.65\,b_{0.25,0.15}$ and $\beta_2=0.9\,b_{0.81,0.095}$:

$$
H(x)=\big(\beta_1(x)-0.26\big)_+ + 0.72\,\beta_2(x),
$$
$$
g(t)=(2e^{-t})^{\kappa}\Big[\big(0.65\,b_{0.25,0.15}(2e^{-t}-1)-0.26\big)_+ + 0.648\,b_{0.81,0.095}(2e^{-t}-1)\Big].
$$

**Exact error, proved for every $x$.** Since $u_g(x,T)=H(x)$ exactly, we have $|H-D|\le0.26$ on the first peak and $|H-D|\le0.252$ on the second. So $E(g)=0.26\le 0.29$.

**Checks of the restrictions.**

| Check | Analytic result | Notebook estimate |
|---|---|---|
| $g(0),\ g(T)$ | both 0 (outside the supports) | 0, 0 |
| $\min g$ | 0 | 0 |
| $\max g$ | $\le 0.909 < 1$ | 0.8849 |
| $\int_0^T g$ | $\le 0.0862 \le 0.093$ | 0.0837 |
| $E(g)$ | exactly 0.26 | 0.260000 |
| Smoothness | continuous, piecewise $C^1$, two corners at $t=0.4045,\ 0.5401$ | — |

---

## 5. Part C: The numerical method

### C1. Update rule and positivity

The supplied solver uses a forward difference in time and a backward (upwind) difference in space:

$$
U_j^{n+1}=(1-\kappa\Delta t-\lambda_j)\,U_j^n+\lambda_j\,U_{j-1}^n,\qquad \lambda_j=\frac{(1+x_j)\Delta t}{h}.
$$

The time step is chosen as $M=\lceil T(\kappa+2N)/0.8\rceil$, which gives $\kappa\Delta t+\lambda_j\le0.8$. Consequences:

- Both coefficients are nonnegative (the first is at least 0.2).
- By induction, $U_j^n\ge0$ whenever $g\ge0$.
- The coefficients sum to at most 1, so $U\le\max g$: there is no overshoot.

### C2. Results and the effect of refinement

By the triangle inequality, $E_{\mathrm{comp},N}\le E(g)+\delta_N$.

| Target | $N$ | $M$ | $E(g)$ | $\delta_N$ | $E_{\mathrm{comp},N}$ | Computed / exact peak at $x^*$ | Nodes on peak 2 |
|---|---|---|---|---|---|---|---|
| original ($w=0.095$) | 100 | 174 | 0.2600 | 0.3759 | 0.6277 | 0.275 / 0.648 | 19 |
| | 200 | 348 | 0.2600 | 0.2893 | 0.5412 | 0.360 / 0.648 | 38 |
| | 400 | 694 | 0.2600 | 0.1988 | 0.4507 | 0.450 / 0.648 | 76 |
| narrow ($w_n=0.025$) | 100 | 174 | 0.2600 | 0.5711 | 0.8231 | 0.078 / 0.648 | 5 |
| | 400 | 694 | 0.2600 | 0.4961 | 0.7481 | 0.153 / 0.648 | 20 |

**Observations:**

- **Refinement helps, but slowly.** Each halving of $h$ reduces $\delta_N$ by about 23% and then 31%. Extra runs give $0.122$ at $N=800$ and $0.068$ at $N=1600$, so the first-order rate appears only on fine grids.
- **Most of the error comes from smearing of peak 2.** The upwind averaging acts like numerical diffusion, $\approx\frac h2(1+x)(1-\lambda_j)u_{xx}$. Peak 2 enters the domain early as a narrow pulse and gets averaged during its whole trip of length $\log 1.81 = 0.59$.
- **The bound is strict.** $E(g)$ is reached on peak 1, while $\delta_N$ peaks near $x^*$. So $E_{\mathrm{comp},N}\approx 0.252+\delta_N<E(g)+\delta_N$.
- **The design affects the numerics.** Keeping peak 2 wide reduced $\delta_{400}$ from 0.255 (candidate ii) to 0.199.

---

## 6. Part D: Narrower target

### D1. Prediction and construction

**Prediction.** The peak height is unchanged, so the cap problem stays the same. The boundary pulse becomes about 4× shorter, so the budget is no longer the binding restriction. Numerically, the narrow peak should be smeared much more, giving a larger $\delta_N$.

**Construction.** Since the budget is relaxed, I scale both peaks without any clipping. This gives a boundary function that is $C^1$ everywhere:

$$
g_n(t)=(2e^{-t})^{\kappa}\big[0.39\,b_{0.25,0.15}(2e^{-t}-1)+0.648\,b_{0.81,0.025}(2e^{-t}-1)\big].
$$

**Checks:**

| Quantity | Result |
|---|---|
| $E(g_n)$ | exactly 0.26 |
| $\max g_n$ | $\le 0.891$ (estimate 0.8848) |
| $\int g_n$ | $\le 0.0682$ (estimate 0.0649) |
| $g_n(0),\ g_n(T)$ | 0, 0 |
| Smoothness | $C^1$, no corners |

### D2. Comparison with the original target

- **Integral:** $\int g$ drops from 0.0837 to 0.0649. **Exact error:** $E=0.26$ for both targets, as predicted.
- **Numerical error:** $\delta_N$ is much larger: 0.571 vs 0.376 at $N=100$, and 0.496 vs 0.199 at $N=400$. Refining from $N=100$ to $400$ improves it by only 13%, against 47% for the original target. Even $N=1600$ leaves $\delta_N=0.363$.
- **Explanation:** the narrow peak contains only 5 nodes at $N=100$ and 20 at $N=400$, compared with 19 and 76 for the original. The upwind averaging flattens it to 0.078 and 0.153, against the exact value 0.648.

---

## 7. Conclusions

- The method of characteristics turns the design problem into choosing the final profile: $u_g(x,T)=H(x)$ exactly when $g(s)=(2e^{-s})^{\kappa}H(2e^{-s}-1)$.
- Damping together with the cap $g\le1$ forces $E\ge0.168$, and the budget also rules out exact reproduction.
- My per-peak design meets every restriction with **$E=0.26\le0.29$**, proved for all $x$, for both targets.
- The upwind scheme is positive and stable but diffusive. It underestimates narrow peaks, and that numerical error, not the design, is why $E_{\mathrm{comp},N}>\varepsilon$.

---

## 8. Figures

| Figure | Content |
|---|---|
| `partA_characteristics.png` | Characteristic curves, the two regions, and the $t_0=0.45$ line |
| `figB_design.png` | $g_D$ vs. the chosen $g$ and $g_n$; the exact final profile inside the band $D\pm\varepsilon$ |
| `original_right.png`, `narrow_right.png` | Notebook output: target, exact and computed solutions |
| `figC_refinement.png` | $\delta_N$ vs. $N$ (log–log) and a zoom on the second peak at $N=400$ |

---

## 9. How to reproduce the results

1. Open `MSO203M_MiniProject1_250120_Amritanshu_Kumar.ipynb` in [Google Colab](https://colab.research.google.com).
2. Run Step 1 to load the tools.
3. In Step 2, set `roll_number = "250120"`.
4. In Step 3, enter the original-target formula:
   ```
   (2*exp(-t))^kappa*(pos(0.65*bump(2*exp(-t)-1,0.25,0.15)-0.26)+0.72*0.9*bump(2*exp(-t)-1,xstar,w))
   ```
5. In Step 4, enter the narrower-target formula:
   ```
   (2*exp(-t))^kappa*(0.6*0.65*bump(2*exp(-t)-1,0.25,0.15)+0.72*0.9*bump(2*exp(-t)-1,xstar,w))
   ```
6. Run Step 5 to export the results ZIP.
7. Paste `colab_figures_cell.py` into a new cell and run it to regenerate all the report figures.

---

## 10. Repository structure

```
├── README.md
├── Mini_Project_1_Report_Amritanshu_250120.pdf          # 4-page report (Parts A–D)
├── MSO203M_MiniProject1_250120_Amritanshu_Kumar.ipynb   # notebook with all runs and outputs
├── MSO203M_250120_results.zip                           # results exported by notebook Step 5
├── colab_figures_cell.py                                # code to regenerate the figures
├── report/
│   └── main.tex                                         # LaTeX source of the report
└── figures/                                             # PNG figures used in the report
```

---

## 11. Tools used

Python (NumPy, Matplotlib), Google Colab / Jupyter, LaTeX (Overleaf)

---

## 12. Acknowledgements

- The project brief, the reading handout and the notebook's solver and checking code were provided by the MSO203M course instructor, IIT Kanpur.
- I used Claude (Anthropic) as an assistant for the derivations, for comparing candidate designs and for drafting the report, as the course allowed. I checked all results by hand and by running the notebook myself.
