# Note: why ConstantKernel / bounds / noise / RBF are all "dials"

Context: slide 50-51, Example 2. Line in question:

```python
kernel = ConstantKernel(2) * RBF(10, length_scale_bounds=(1, 1e2))  # "fixed length scale"
```

All numbers below are measured on `example2Slide50Data.txt` with `alpha=0.2`,
`normalize_y=True`. Reproduce with `/tmp/dials.py`.

---

## 1. The numbers you type are seeds, not settings

That line fits to:

```
1.12**2 * RBF(length_scale=1.29)
```

**Neither number survived.** The `2` became 1.25, the `10` became 1.29.
`gp.fit()` runs L-BFGS on the log marginal likelihood and overwrites both.
So the comment "fixed length scale" is describing the opposite of what happens.

Starting value is nearly irrelevant when bounds are open — l0 = 0.5, 10, and 90
all converge to the same 1.29:

| start l | fitted l | LML |
|---|---|---|
| 0.5 | 1.29 | -20.76 |
| 10 | 1.29 | -20.76 |
| 90 | 1.29 | -20.76 |

To actually fix a value you must say so:

```python
RBF(10, length_scale_bounds="fixed")   # truly fixed at 10
```

## 2. The real dial is the BOUNDS, not the value

Bounds are where the optimizer is allowed to look. Move the floor up and the fit
gets clamped at it:

| bounds | fitted l | LML |
|---|---|---|
| (1, 100) | 1.29 | **-20.76** |
| (3, 100) | 3.00 | -27.42 |
| (10, 100) | 10.0 | -43.34 |

Your `(1, 1e2)` happens to be safe — the true optimum 1.29 sits just above the
floor of 1. Had you written `(10, 1e2)` you'd have silently lost 22 LML.

sklearn *does* warn on this:

> ConvergenceWarning: The optimal value found for dimension 0 of parameter
> k2__length_scale is close to the specified lower bound 10.0.

but it's one line in a wall of output — the slide-51 sweep emitted a dozen of
these and they were easy to scroll straight past.

## 3. The dials are NOT independent — they compensate

This is the real reason "they're all dials" is the wrong mental model. Watch the
signal variance inflate as the length scale is forced up:

| forced l | ConstantKernel becomes |
|---|---|
| 1.29 | 1.12² |
| 3 | 2.54² |
| 10 | 7.60² |

Squeeze one, another bulges. A long length scale makes the GP too smooth to
track the bump at x=5, so the optimizer buys back flexibility by cranking the
amplitude. They are coupled parameters of one prior, not four knobs on a panel.

## 4. What each one actually means

- **`RBF(length_scale)`** — correlation distance. How far apart two x's must be
  before the GP treats them as unrelated. *The dominant dial.* Too short →
  chases noise; too long → flattens through the bump.
- **`ConstantKernel(v)`** — signal variance (amplitude). How far the function is
  allowed to roam from its mean. Weak dial when free; one-sided when fixed
  (too small badly damps the fit, too large is harmless — see slide 51 Q4).
- **`alpha`** — noise variance on the diagonal. Does double duty: it's your
  measurement-noise model *and* numerical jitter keeping the covariance matrix
  invertible. Trades honoring the data against smoothing through it.
- **`*_bounds`** — the optimizer's search box, and a crude prior. Usually the
  thing you actually want to set.

## 5. Absolute scale matters, not just the signal-to-noise ratio

Tempting shortcut: "only signal/noise ratio matters." **Measured, it doesn't hold.**
Three kernels at ratio 5, length scale fixed at 1.289:

| const | alpha | ratio | total prior var | LML |
|---|---|---|---|---|
| 1 | 0.2 | 5 | 1.2 | **-20.83** |
| 2 | 0.4 | 5 | 2.4 | -23.85 |
| 10 | 2.0 | 5 | 12.0 | -37.88 |

Same ratio, 17 LML apart. Because `normalize_y=True` standardizes y to unit
variance, so **signal + noise should total roughly 1** — the data's variance.
`const=1, alpha=0.2` wins because 1.2 ≈ 1. The ratio sets how much structure is
called signal vs noise; the sum has to match the data's actual spread.

---

**Takeaway:** the value you type is a guess, the bounds are the setting, and the
parameters trade against each other. Check `gp.kernel_` after every fit — it is
the only statement of what the model actually used.
