# Sundials

Notes on Sundials library

## CVode

We consider the general $N$-dimension IVP problem

$$
    \dot{y} = f(t, y) \qquad y(t_0) = y_0
$$
where $y \in \R^N$ and $f: \R \times \R^N \to \R^N$.

The methods used in CVode are variable order, variable multi-step methods of the form

$$
    \sum_{i=0}^{K_1} \alpha_{n, i} y^{n - i} + h_n \sum_{i=0}^{K_2} \beta_{n, i} \dot{y}^{n-i} = 0
$$
where $y_n$ is the numerical approximation of $y(t_n)$ and $h_n = t_n - t_{n-1}$ is the step size.

Two types of methods are available in CVode:

1. Adams-Moulton: $K_1 = 1, K_2 = q - 1, \: q \in \intint{1}{12}$
2. Backward Differentiation Methods: $K_1 = q, K_2 = 0, \: q \in \intint{1}{5}$

The coefficients are determined by the method type, the order, the recent history of step sizes with the normalization coefficient $\alpha_{n, 0} = - 1$.

For either choice of formula, a non-linear system must be solved of the form

$$
    F(y^n) = y^n - h_n \beta_{n, 0} f(t_n, y^n) - a_n = 0
$$

where $a_n = \sum_{i>0} [\alpha_{n, i} y^{n - i} + h_n \beta_{n, i} \dot{y}^{n-i}]$.

Resolution is done using a Newton method:

$$
    M[y^{n(m + 1)} - y^{n(m)}] = -F(y^{n(m)})
$$

where $M = I - \gamma J, J = \pdv{f}{y}, \gamma = h_n \beta_{n, 0}$.
