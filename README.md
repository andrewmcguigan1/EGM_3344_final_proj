# EGM 3344 Final Project (Fall 2025)
## Polynomial regression of airfoil lift data + numerical stability of least-squares solvers

This repository contains a MATLAB implementation of a regression problem (fit lift coefficient vs.
angle of attack) and a numerical-analysis comparison of several least-squares solution approaches.

The core theme is **how algorithm choice and problem conditioning affect accuracy and stability**
as the regression model becomes more complex (higher polynomial degree).

---

## 1) Motivation

In aerodynamic design and analysis, it is common to model a relationship such as lift coefficient
$C_L$ versus angle of attack $\alpha$ using data from experiments or simulations. A simple linear
model can be useful, but real data often shows curvature and nonlinearity. Polynomial regression is
an easy way to increase model flexibility, but (unfortunately for everyone involved) higher-degree
polynomials can become **numerically unstable**.

This project uses a small dataset of $(\alpha, C_L)$ measurements for NACA 24xx airfoils and explores:

1. How well different polynomial degrees fit the data.
2. How to solve the resulting least-squares problem using different numerical methods.
3. How solution quality degrades as the polynomial degree increases and the design matrix becomes
   ill-conditioned.

---

## 2) Repository contents

| File | Purpose |
|---|---|
| `Project.m` | Main driver script: builds regression models, runs comparisons, prints tables, and generates plots. |
| `givens.m` | Custom QR factorization for least squares using **Givens rotations** (returns $R$ and $Q^T b$). |
| `qrfact.m` | QR factorization using **Householder reflections** (provided “from textbook” style). |
| `final_proj_instructions.pdf` | Course project description / requirements handout. |

---

## 3) Mathematical formulation

### 3.1 Data and regression model

Let there be $m$ measured data points $(\alpha_i, y_i)$ where:

* $\alpha_i$ is the angle of attack (degrees)
* $y_i$ is the measured lift coefficient $C_L$

For a polynomial model of degree $n$, we use

$$
\hat{y}(\alpha) = p_n(\alpha) = \beta_0 + \beta_1\alpha + \beta_2\alpha^2 + \cdots + \beta_n\alpha^n.
$$

Define the (Vandermonde-like) design matrix $X \in \mathbb{R}^{m\times (n+1)}$ as

$$
X_{i,k} = \alpha_i^{k-1}, \quad k = 1,\dots,n+1.
$$

Then the model predictions are

$$
\hat{y} = X\beta,
$$

where $\beta = [\beta_0,\beta_1,\dots,\beta_n]^T$.

### 3.2 Least-squares problem

Because the system is typically overdetermined ($m > n+1$), we solve the **linear least-squares**
problem:

$$
\min_{\beta}\ \|X\beta - y\|_2.
$$

The numerical question is not “can we solve it?” but rather “can we solve it **without making the
computer cry** when $X$ becomes ill-conditioned?”

---

## 4) Algorithms compared

### 4.1 QR factorization (reference)

If $X = QR$ where $Q$ has orthonormal columns and $R$ is upper-triangular, then the least-squares
solution satisfies

$$
\beta = R^{-1}Q^T y.
$$

This project uses:

* MATLAB’s built-in `qr(X,0)` as a reference baseline.
* A Householder-based QR (`qrfact.m`).
* A Givens-rotation QR specialized to compute $R$ and $Q^T y$ directly (`givens.m`).

Why QR? Orthogonal transformations preserve the 2-norm:

$$
\|Qz\|_2 = \|z\|_2,
$$

which is a big reason QR methods are typically stable for least squares.

### 4.2 Givens rotations (implemented in `givens.m`)

Givens rotations eliminate individual subdiagonal entries with a 2D rotation. For entries
$a = R_{j,j}$ and $g = R_{i,j}$, choose $c$ and $s$ so that:

$$ \begin{bmatrix} c & -s \\\\ s & c \end{bmatrix} \begin{bmatrix} a \\\\ g \end{bmatrix} = \begin{bmatrix} r \\\\ 0 \end{bmatrix}, \quad r = \sqrt{a^2 + g^2},\quad c = \frac{a}{r},\quad s = -\frac{g}{r}. $$

Applying these rotations across columns transforms $X$ into an upper-triangular $R$. Applying the
same rotations to the right-hand side $y$ yields $Q^T y$ without explicitly forming $Q$.

### 4.3 Householder reflections (implemented in `qrfact.m`)

Householder QR eliminates entire subcolumns at once using reflectors

$$
H = I - 2vv^T,
$$

with a carefully chosen vector $v$ to zero out entries below the diagonal.

### 4.4 Normal equations + Cholesky

Another approach is to form the normal equations:

$$
X^T X\,\beta = X^T y.
$$

If $X^T X$ is symmetric positive definite, we can solve via Cholesky

$$
X^T X = R^T R.
$$

However, conditioning gets worse because (roughly)

$$
\kappa_2(X^T X) = \kappa_2(X)^2,
$$

which is why Cholesky on the normal equations often breaks down earlier than QR for poorly
conditioned problems.

---

## 5) Metrics used in the comparisons

The script evaluates methods using:

### 5.1 Fit quality

* **RMSE**:  
$\mathrm{RMSE} = \sqrt{\frac{1}{m}\sum_{i=1}^{m}(y_i - \hat{y}_i)^2}.$

* **Coefficient of determination** ($R^2$):
  $$
  R^2 = 1 - \frac{\sum_{i=1}^{m}(y_i-\hat{y}_i)^2}{\sum_{i=1}^{m}(y_i-\bar{y})^2},
  \quad \bar{y} = \frac{1}{m}\sum_{i=1}^{m}y_i.
  $$

### 5.2 Conditioning / stability indicators

* **Condition number** (2-norm):
  $$
  \kappa_2(X) = \|X\|_2\,\|X^{-1}\|_2.
  $$

* **Relative coefficient error** vs MATLAB QR:
  $$
  \frac{\|\beta_{\text{method}} - \beta_{\text{QR}}\|_2}{\|\beta_{\text{QR}}\|_2}.
  $$

### 5.3 Runtime and flop estimates

`Project.m` measures runtime using `tic/toc` and uses a rough flop estimate used in the script:

$$
\text{flops} \approx 6mc - 3c(c+1),
$$

where $m$ is the number of data points and $c$ is the number of columns in $X$ (for a degree-$n$
polynomial, $c=n+1$). This is an *approximation* meant for comparison, not a perfect instruction-by-
instruction count.

---

## 6) What `Project.m` does (high-level)

### Part A: Givens rotation regression (degrees 1 through 6)

1. Fits a linear model ($n=1$) using the custom Givens QR.
2. Fits polynomial models for degrees $n=2$ to $6$ (Vandermonde matrix).
3. For each fit, computes RMSE, $R^2$, condition number, runtime, and flop estimate.
4. Produces comparison plots for each polynomial degree.

### Part B: Leave-one-out cross validation (LOOCV) for degrees 1 through 24

For each degree $n$ and each data point $j$:

1. Remove point $j$.
2. Fit coefficients $\beta^{(-j)}$ using Givens QR.
3. Predict the left-out point:
   $$
   \hat{y}_j = x_j^T\beta^{(-j)}.
   $$

It then reports the LOOCV RMSE vs polynomial degree and computes variance of fitted coefficients
across the LOOCV runs (a rough indicator of coefficient stability).

### Part C: Stability comparison vs degree (QR vs Householder vs Cholesky vs Givens)

For degrees $n=2$ to $25$:

1. Compute $\kappa_2(X)$.
2. Solve for $\beta$ using:
   * MATLAB `qr`
   * Householder (`qrfact`)
   * Cholesky on normal equations
   * Givens (`givens`)
3. Compare RMSE and relative coefficient error vs the MATLAB QR result.

---

## 7) Requirements and how to run

### Requirements

* MATLAB (recommended). The script uses:
  * `table` for formatted outputs
  * `cond`, `qr`, `chol`
  * basic plotting (`plot`, `semilogy`)

No special toolboxes are required beyond core MATLAB functionality.

### Running

1. Open MATLAB.
2. Set your current folder to the repository directory.
3. Run:

```matlab
Project
```

### Selecting a dataset

In `Project.m`, choose which airfoil data vector to fit by editing:

```matlab
y = CL_2412;  % or CL_2414, CL_2415
```

---

## 8) Expected outputs

When you run `Project.m`, you should see:

* Multiple figures comparing the raw data, the linear fit, and polynomial fits (degrees 2 through 6).
* Printed tables in the MATLAB console for:
  * Linear fit summary
  * Polynomial fit summaries
  * LOOCV RMSE vs degree
  * Method comparison table (conditioning, RMSE, relative coefficient errors)
* Two semilog plots:
  * $\kappa_2(X)$ vs polynomial degree
  * relative coefficient error vs polynomial degree (method breakdown)

---

## 9) Practical insights (what you’re supposed to learn here)

1. **Higher-degree polynomials are not “free accuracy.”** They may reduce training error while
   increasing generalization error (overfitting). LOOCV helps reveal this tradeoff.

2. **Vandermonde matrices become ill-conditioned quickly.** As degree increases, the columns
   (powers of $\alpha$) become increasingly correlated, causing $\kappa_2(X)$ to grow and making the
   computed coefficients sensitive to rounding error.

3. **Normal equations can be numerically fragile.** Since $\kappa(X^T X) \approx \kappa(X)^2$,
   Cholesky on the normal equations tends to deviate from the QR-based reference sooner.

4. **Orthogonal methods (Householder / Givens) are typically more stable.** Both approaches are
   based on orthogonal transformations, which generally provide better numerical behavior for least
   squares than explicitly forming $X^T X$.

5. **If you care about robustness, scale the problem.** A common improvement is to scale/shift
   inputs (e.g., map $\alpha$ to $[-1,1]$) or use orthogonal polynomial bases (Chebyshev/Legendre)
   instead of raw monomials $\alpha^k$.

---

## 10) Notes / reproducibility

* `tic/toc` timings depend on machine load. For more reliable benchmarking, run multiple trials and
  average.
* The flop estimate is a coarse comparison tool, not a cycle-accurate performance model.
* If Cholesky fails for some degrees, `Project.m` catches the failure and records `NaN`.
