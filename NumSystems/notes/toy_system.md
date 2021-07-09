# Toy system for linear systems

Let us consider the simple 1D Dirichlet Poisson problem:

$$
\begin{align}
\nabla^2 \phi &= \frac{\partial^2 \phi}{\partial x^2} = -R \ \ \text{in} \ \dot{\Omega} \\
\phi &= \phi_D \ \ \text{on} \ \partial \Omega
\end{align}
$$

The system yields for $n$ nodes:

$$
\begin{bmatrix}
1 & 0 & \ldots & \ldots & \ldots & 0 \\
1 & -2 & 1 & \ldots & \ldots & 0 \\
0 & 1 & -2 & 1 & \ldots & 0 \\
0 & 0 & \ddots & \ddots & \ldots & 0 \\
0 & 0 & 0 & 1 & -2 & 1 \\
0 & 0 & 0 & \ldots & 0 & 1
\end{bmatrix}
\begin{bmatrix}
\phi_1 \\ \vdots \\ \\ \\ \vdots\\ \phi_n
\end{bmatrix} =
\begin{bmatrix}
\phi_l \\ - R_2 \Delta x^2 \\ \vdots \\ \vdots \\ - R_{n-1} \Delta x^2 \\ \phi_r
\end{bmatrix}
$$

The inverse of the problem matrix is

$$
\begin{bmatrix}
1 & 0 & \ldots & \ldots & 0 & 0 \\
\frac{n-2}{n-1} & -2 & 1 & \ldots & \ldots & \frac{1}{n-1} \\
\frac{n-3}{n-1} & 1 & -2 & 1 & \ldots & \frac{2}{n-1} \\
\vdots & 0 & \ddots & \ddots & \ldots & \vdots \\
\frac{1}{n-1} & 0 & 0 & 1 & -2 & \frac{n-2}{n-1} \\
0 & 0 & \ldots & \ldots & 0 & 1
\end{bmatrix}
$$

The analytical solution of the problem is:

$$
\begin{align}
\phi(x) &= \sum_{n=1}^{+\infty} \phi_{n} \sin(n \pi x / L_x) \\
\phi_n &= \left(\frac{L_x}{n\pi}\right)^2 \frac{2}{L_x}\int_0^{L_x} \sin(n \pi x' / L_x) R(x') \mrm{d}x'
\end{align}
$$