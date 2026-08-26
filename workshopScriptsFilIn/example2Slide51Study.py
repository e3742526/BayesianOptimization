"""Slide 51 'Lets Try It' -- parameter study for Example 2.

Answers the four questions on slide 51:
  1. What noise variance gives the best fit?
  2. What length scale gives the best fit with alpha = 0.2?
  3. What length scale gives the worst fit with alpha = 0.2?
  4. How sensitive is the output to the signal noise?

Fit quality is scored two ways:
  - LML     : log marginal likelihood (the GP's own model-selection criterion)
  - LOO-CV  : leave-one-out RMSE (intuitive 'does it predict held-out points')
Both are reported because they can disagree, and the disagreement is the lesson.
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, ConstantKernel, WhiteKernel

OUT = Path(__file__).resolve().parents[1] / "output"
OUT.mkdir(exist_ok=True)

data = np.loadtxt(Path(__file__).with_name("example2Slide50Data.txt"), delimiter=",")
X = data[:, 0].reshape(-1, 1)
y = data[:, 1]


def fit(length_scale, alpha, const=1.0, free_const=True, free_ls=False):
    """Fit a GP. Hyperparameters are fixed unless explicitly freed."""
    c = ConstantKernel(
        const, constant_value_bounds=(1e-3, 1e3) if free_const else "fixed"
    )
    r = RBF(length_scale, length_scale_bounds=(1e-2, 1e3) if free_ls else "fixed")
    kernel = c * r
    gp = GaussianProcessRegressor(
        kernel=kernel,
        alpha=alpha,
        normalize_y=True,
        optimizer="fmin_l_bfgs_b" if (free_const or free_ls) else None,
        n_restarts_optimizer=2 if (free_const or free_ls) else 0,
    )
    gp.fit(X, y)
    return gp


def loo_rmse(length_scale, alpha, const):
    """Leave-one-out RMSE at FIXED hyperparameters (no refitting per fold)."""
    errs = []
    for i in range(len(X)):
        m = np.ones(len(X), dtype=bool)
        m[i] = False
        g = GaussianProcessRegressor(
            kernel=ConstantKernel(const, constant_value_bounds="fixed")
            * RBF(length_scale, length_scale_bounds="fixed"),
            alpha=alpha,
            normalize_y=True,
            optimizer=None,
        )
        g.fit(X[m], y[m])
        errs.append(y[i] - g.predict(X[i].reshape(1, -1))[0])
    return float(np.sqrt(np.mean(np.square(errs))))


print("=" * 68)
print("Q1: What noise variance gives the best fit?")
print("=" * 68)

# Let the GP learn the noise itself via a WhiteKernel -- the direct answer.
gp_white = GaussianProcessRegressor(
    kernel=ConstantKernel(1.0) * RBF(1.0) + WhiteKernel(0.1),
    normalize_y=True,
    n_restarts_optimizer=5,
)
gp_white.fit(X, y)
learned_noise = gp_white.kernel_.k2.noise_level
print(f"  WhiteKernel-learned noise variance : {learned_noise:.4f}")
print(f"  (slide uses alpha = 0.2)           : ratio {0.2 / learned_noise:.1f}x too big")
print(f"  learned kernel: {gp_white.kernel_}")

alphas = np.logspace(-3, 0.3, 40)
a_lml, a_loo = [], []
for a in alphas:
    g = fit(1.0, a, free_const=True, free_ls=True)
    a_lml.append(g.log_marginal_likelihood_value_)
    ls = g.kernel_.k2.length_scale
    cc = g.kernel_.k1.constant_value
    a_loo.append(loo_rmse(ls, a, cc))
a_lml, a_loo = np.array(a_lml), np.array(a_loo)
best_a_lml = alphas[int(np.argmax(a_lml))]
best_a_loo = alphas[int(np.argmin(a_loo))]
print(f"  best alpha by LML                  : {best_a_lml:.4f}")
print(f"  best alpha by LOO-CV RMSE          : {best_a_loo:.4f}")

print()
print("=" * 68)
print("Q2/Q3: Best and worst length scale at alpha = 0.2")
print("=" * 68)

ALPHA = 0.2
lss = np.logspace(-1, 1.4, 60)
l_lml, l_loo, l_const = [], [], []
for ls in lss:
    g = fit(ls, ALPHA, free_const=True, free_ls=False)
    cc = g.kernel_.k1.constant_value
    l_lml.append(g.log_marginal_likelihood_value_)
    l_const.append(cc)
    l_loo.append(loo_rmse(ls, ALPHA, cc))
l_lml, l_loo = np.array(l_lml), np.array(l_loo)

best_l_lml = lss[int(np.argmax(l_lml))]
worst_l_lml = lss[int(np.argmin(l_lml))]
best_l_loo = lss[int(np.argmin(l_loo))]
worst_l_loo = lss[int(np.argmax(l_loo))]

print(f"  BEST  length scale by LML     : {best_l_lml:.3f}  (LML {l_lml.max():.2f})")
print(f"  BEST  length scale by LOO-CV  : {best_l_loo:.3f}  (RMSE {l_loo.min():.4f})")
print(f"  gp.fit() MLE from the slide 50 script : 1.289")
print(f"  WORST length scale by LML     : {worst_l_lml:.3f}  (LML {l_lml.min():.2f})")
print(f"  WORST length scale by LOO-CV  : {worst_l_loo:.3f}  (RMSE {l_loo.max():.4f})")
print(f"  grid spanned {lss[0]:.2f} .. {lss[-1]:.2f}")

print()
print("=" * 68)
print("Q4: Sensitivity to the signal variance (ConstantKernel)")
print("=" * 68)

consts = np.logspace(-2, 2, 25)
c_lml, c_loo, c_peak = [], [], []
Xs = np.linspace(0, 10, 200).reshape(-1, 1)
for c in consts:
    g = fit(1.289, ALPHA, const=c, free_const=False, free_ls=False)
    c_lml.append(g.log_marginal_likelihood_value_)
    c_loo.append(loo_rmse(1.289, ALPHA, c))
    c_peak.append(float(g.predict(Xs).max()))
c_lml, c_loo, c_peak = np.array(c_lml), np.array(c_loo), np.array(c_peak)
print(f"  signal variance {consts[0]:.3g} .. {consts[-1]:.3g}")
print(f"  LOO RMSE range  : {c_loo.min():.4f} .. {c_loo.max():.4f}")
print(f"  peak height     : {c_peak.min():.3f} .. {c_peak.max():.3f}  (data peak {y.max():.3f})")
plateau = consts[c_loo < c_loo.min() * 1.05]
if len(plateau):
    print(f"  within 5% of best LOO for signal var >= {plateau.min():.3g}")

# ----------------------------------------------------------------- figures
fig, ax = plt.subplots(2, 2, figsize=(15, 10))

ax[0, 0].semilogx(lss, l_lml, "b-", lw=2)
ax[0, 0].axvline(best_l_lml, color="g", ls="--", label=f"best {best_l_lml:.2f}")
ax[0, 0].axvline(1.289, color="k", ls=":", label="gp.fit() 1.29")
ax[0, 0].set_xlabel("length scale")
ax[0, 0].set_ylabel("log marginal likelihood")
ax[0, 0].set_title("Q2/Q3: LML vs length scale (alpha=0.2)")
ax[0, 0].legend()

ax[0, 1].semilogx(lss, l_loo, "r-", lw=2)
ax[0, 1].axvline(best_l_loo, color="g", ls="--", label=f"best {best_l_loo:.2f}")
ax[0, 1].set_xlabel("length scale")
ax[0, 1].set_ylabel("LOO-CV RMSE")
ax[0, 1].set_title("Q2/Q3: LOO-CV vs length scale")
ax[0, 1].legend()

ax[1, 0].semilogx(alphas, a_lml, "b-", lw=2)
ax[1, 0].axvline(best_a_lml, color="g", ls="--", label=f"best {best_a_lml:.3f}")
ax[1, 0].axvline(0.2, color="k", ls=":", label="slide 0.2")
ax[1, 0].axvline(learned_noise, color="m", ls="-.", label=f"white {learned_noise:.3f}")
ax[1, 0].set_xlabel("alpha (noise variance)")
ax[1, 0].set_ylabel("log marginal likelihood")
ax[1, 0].set_title("Q1: LML vs noise variance")
ax[1, 0].legend()

ax[1, 1].semilogx(consts, c_loo, "r-", lw=2)
ax[1, 1].set_xlabel("signal variance (ConstantKernel)")
ax[1, 1].set_ylabel("LOO-CV RMSE")
ax[1, 1].set_title("Q4: sensitivity to signal variance")

plt.tight_layout()
plt.savefig(OUT / "example2-slide51-sweeps.png", dpi=110)
print(f"\nwrote {OUT / 'example2-slide51-sweeps.png'}")

# best vs worst fit, side by side
fig2, axes = plt.subplots(1, 3, figsize=(18, 5))
for a, ls, title in [
    (axes[0], 0.15, "too short -> chases noise"),
    (axes[1], best_l_lml, f"best  l={best_l_lml:.2f}"),
    (axes[2], 15.0, "too long -> flat, misses bump"),
]:
    g = fit(ls, ALPHA, free_const=True, free_ls=False)
    m, s = g.predict(Xs, return_std=True)
    a.plot(X, y, "rs", markersize=9)
    a.fill_between(Xs.ravel(), m - s, m + s, alpha=0.3)
    a.plot(Xs, m, "b-", lw=2)
    a.set_title(title)
    a.set_xlabel("Parameter")
    a.set_ylim(-0.8, 1.4)
axes[0].set_ylabel("Blast Resistance")
plt.tight_layout()
plt.savefig(OUT / "example2-slide51-bestworst.png", dpi=110)
print(f"wrote {OUT / 'example2-slide51-bestworst.png'}")
