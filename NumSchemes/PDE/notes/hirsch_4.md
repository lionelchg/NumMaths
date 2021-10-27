# Numerical Computation of Internal and External Flows notes

## Chapter 8: General properties and High-Resolution Numerical Schemes

### Finite volume formulation of schemes and limiters

#### Cell-face interpolation

Considering interpolation on the cell values in 1D

$$
u_{i+1/2} = u_i + \alpha_i(u_{i + 1} - u_i) + \beta_i(u_i - u_{i - 1})
$$

Expanding in Taylor expansion around $i + 1/2$

$$
u_{i+1/2} = u_i + \frac{\Delta x_{i+1}}{2} (u_x)_i + \frac{\Delta x_{i+1}^2}{8}(u_{xx})_i + \mc{O}(\Delta x^3)
$$

Taylor expansion of the right hand side yields

$$
u_i + \alpha_i (\Delta x_{i+1} (u_x)_i + \frac{\Delta x_{i+1}^2}{2}(u_{xx})_i) + \beta_i (\Delta x_i (u_x)_i - \frac{\Delta x_i^2}{2}(u_{xx})_i)
$$

Two conditions for third order accuracy are

$$
\begin{aligned}
\alpha_i \Delta x_{i+1} + \beta_i \Delta x_i = \frac{\Delta x_{i+1}}{2}\\
\alpha_i \Delta x_{i+1}^2 - \beta_i \Delta x_i^2 = \frac{\Delta x_{i+1}^2}{4}
\end{aligned}
$$

So that for regular meshes and second order accuracy

$$
\alpha + \beta = \frac{1}{2}
$$

#### Numerical flux

We consider the one-dimensional conservation law

$$
\pdv{u}{t} + \pdv{f}{x} = 0
$$

Discretization of the conservation law in finite volume formulation yields

$$
\frac{\mrm{d}u_i}{\mrm{d}t} = -\frac{1}{\Delta x_i}(f^*_{i+1/2} - f^*_{i-1/2})
$$

where the numerical flux at cell $i + 1/2$ depends on the surrounding points

$$
f^*_{i+1/2}((u_k)_{k \in V(i)}).
$$

The consistency condition is

$$
f^*_{i+1/2}((u_k = u)_{k \in V(i)}) = f(u)
$$

For a linear convection problem $f(u) = a u$ so that the interpolated value can be reduced to $u$ in this case.

The general cell-face interpolation is written for second order accuracy and $\beta = (1 - \kappa) / 4$ and is called the $\kappa$-scheme

$$
u_{i+1/2} = u_i + \veps\left[\frac{1 + \kappa}{4}(u_{i+1} - u_i) + \frac{1 - \kappa}{4}(u_i - u_{i-1})\right]
$$

The following table summarizes the remarkable values of this scheme

| Name                 | $\epsilon$ | $\kappa$ |
| -------------------- | ---------- | -------- |
| FOU                  | 0          | -        |
| Second order upwind  | 1          | -1       |
| Central              | 1          | 1        |
| Fromm                | 1          | 0        |
| Quick scheme         | 1          | 1/2      |
| Third order accurate | 1          | 1/3      |

Application of limiter on top of the $\kappa$-scheme is straightforward:

$$
u_{i+1/2} = u_i + \Psi(r_i)\left[\frac{1 + \kappa}{4}(u_{i+1} - u_i) + \frac{1 - \kappa}{4}(u_i - u_{i-1})\right]
$$

However only limiters with $\kappa = -1$ (second order upwind) yield reasonable results. The results with those schemes where $\sig$ is not present at all in the coefficients yield results that are far from the limited finite difference ones.

The best candidate to adapt the scheme into AVIP would be the FV Lax-Wendroff scheme:

$$
u_{i+1/2} = u_i + \frac{1-\sig}{2}(u_{i + 1} - u_i) \Psi(R_i)
$$

where $R_i = (u_i - u_{i - 1}) / (u_{i + 1} - u_i)$