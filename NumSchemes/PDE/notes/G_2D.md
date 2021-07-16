# 2D amplification factors

Notes concerning the derivation of 2D amplification factors. We assume equal spacing $\Delta x$ in $x$ and $y$ direction and we solve

$$
u_t + \vb{a} \cdot \nabla u = 0
$$

with

$$
\vb{a} = a (\cos(\alpha) \vb{e}_x + \sin(\alpha) \vb{e}_y)
$$

Starting from a mode $k$ solution

$$
u(\vb{x}, t) = \hat{V}(k) e^{\I(\vb{k} \cdot \vb{x} - \omega t)}
$$

The exact amplification factor is

$$
\omega = \vb{a} \cdot \vb{k}
$$

## Lax-Wendroff

Starting from a 2D grid with $(i, j)$ indices and making the Taylor expansion in time

$$
u_{ij}^{n+1} = u_{ij}^n + \Delta t (u_t)_{ij}^n + \frac{\Delta t^2}{2} (u_{tt})_{ij}^n
$$

The first time derivative does not raise any issues compared to 1D, the second derivative is a bit more involved:

$$
u_{tt} = \vb{a} \cdot \nabla (\vb{a} \cdot \nabla u) = a_x^2 u_{xx} + a_y^2 u_{yy} + 2 a_x a_y u_{xy}
$$

After some calculus the following relation is derived

$$
\begin{align}
G(\phi_x, \phi_y) &= 1 - \I \sig[\cos(\alpha) \sin(\phi_x) + \sin(\alpha) \sin(\phi_y)] \\
&+ \sig^2[\cos^2\alpha(\cos(\phi_x) - 1) + \sin^2\alpha(\cos(\phi_y) - 1)\\
&+\cos\alpha \sin\alpha(\cos(\phi_x + \phi_y) - \cos(\phi_x - \phi_y))]
\end{align}
$$