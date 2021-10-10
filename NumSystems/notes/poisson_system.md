# Poisson linear system

In these notes we study linear systems associated to the discretization on a mesh of a Poisson problem

$$
\begin{cases}
\nabla^2 \phi = -R \quad &\text{in} \: \dot{\Omega}\\
\phi = \phi_D \quad &\text{on} \: \partial \Omega_D\\
\nabla \phi \cdot \vb{n} = -E_n \quad &\text{on} \: \partial \Omega_N
\end{cases}
$$

The associated linear system will be denoted $Ax = b$.

## Regular mesh

### Cartesian

In 1D cartesian mesh with $x_0 = 0$, $x_N = L$ we have for Dirichlet BCs

$$
A =
\begin{bmatrix}
1 & 0 & \ldots & \ldots & 0 \\
1 & -2 & 1 & \ldots & 0 \\
\vdots & \ddots & \ddots & \ddots & \vdots \\
0 & 0 & 1 & -2 & 1 \\
0 & \ldots & \ldots & 0 & 1 \\
\end{bmatrix}
$$

The matrix can be symmetrized to yield

$$
A =
\begin{bmatrix}
1 & 0 & \ldots & \ldots & 0 \\
0 & -2 & 1 & \ldots & 0 \\
\vdots & \ddots & \ddots & \ddots & \vdots \\
0 & 0 & 1 & -2 & 0 \\
0 & \ldots & \ldots & 0 & 1 \\
\end{bmatrix}
$$

and the removed values will be added to $b$ in this modified version. Let us take the matrix of size $N$ which represents the inner part of the Dirichlet matrices:

$$
A =
\begin{bmatrix}
-2 & 1 &  &  \\
1 & \ddots & \ddots & \\
 & \ddots & \ddots & 1 \\
 &  & 1 & -2 \\
\end{bmatrix}
$$

This matrix has eigenvalues $\lambda_k = 2(1 - \cos(k\pi / (N+1)))$ and eigenvectors $(v_k)_i = \sin(ik\pi / (N+1))$

For Neumann at $x_0$ and Dirichlet at $x_N$ two versions are possible depending on the discretization of the Neumann condition. A ghost-cell approach yields

$$
A =
\begin{bmatrix}
-2 & 2 & \ldots & \ldots & 0 \\
1 & -2 & 1 & \ldots & 0 \\
\vdots & \ddots & \ddots & \ddots & \vdots \\
0 & 0 & 1 & -2 & 1 \\
0 & \ldots & \ldots & 0 & 1 \\
\end{bmatrix}
$$

whereas a first-order discretization gives

$$
A =
\begin{bmatrix}
-1 & 1 & \ldots & \ldots & 0 \\
1 & -2 & 1 & \ldots & 0 \\
\vdots & \ddots & \ddots & \ddots & \vdots \\
0 & 0 & 1 & -2 & 1 \\
0 & \ldots & \ldots & 0 & 1 \\
\end{bmatrix}
$$

### Cylindrical

We consider the 1D Poisson equation in cylindrical coordinates:

$$
\nabla^2 \phi = \frac{1}{r} \frac{\md}{\md r} \left(r \dv{\phi}{r}\right)
$$

The domain goes from $r= 0$ to $r= R$ so that a Neumann boundary condition is mandatory and a Dirichlet BC is applied at the top. For the bottom node special treatment is needed and there is a factor 4.

For interior nodes:

$$
\begin{aligned}
\nabla^2 \phi |_j &= \frac{1}{r_j \Delta r} (r_{j+1/2} \frac{\phi_{j+1} - \phi_j}{\Delta r} - r_{j-1/2} \frac{\phi_j - \phi_{j-1}}{\Delta r}) \\
\nabla^2 \phi |_j &= \frac{1}{\Delta r^2}(\frac{r_{j+1/2} \phi_{j+1} + r_{j-1/2} \phi_{j-1}}{r_j} - 2 \phi_j) \\
\end{aligned}
$$

And for a constant spaced mesh we have therefore

$$
A =
\begin{bmatrix}
-4 & 4 & \ldots & \ldots & 0 \\
1/2 & -2 & 3/2 & \ldots & 0 \\
\vdots & \ddots & \ddots & \ddots & \vdots \\
0 & 0 & (n-3/2)/n & -2 & (n-1/2) / n \\
0 & \ldots & \ldots & 0 & 1 \\
\end{bmatrix}
$$

