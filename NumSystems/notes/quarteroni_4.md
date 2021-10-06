# Chapter 4 Iterative methods for linear systems

We try to find a sequence of vectors that tend to the desired one $Ax = b$

$$
x = \lim_{k \to +\infty} x^{(k)}
$$

Given $x^{(0)}$ an the following iteration procedure is performed:

$$
x^{(k+1)} = B x^{(k)} + f \quad k \geq 0
$$

For consistency we must have $x = Bx + f$ or equivalently $f = (I - B)  A^{-1} b$. This procedure converges if and only if $\rho(B) < 1$. A general method to devise consistent linear method consists in splitting the matrix $A = P - N$ where $P$ is called the preconditioner. From that the iteration is

$$
Px^{(k+1)} = Nx^{(k)} + b \quad k \geq 0
$$

The following relations hold

$$
\begin{aligned}
B &= P^{-1} N \\
x^{(k+1)} &= x^{(k)} + P^{-1} (b - A x^{(k)}) = x^{(k)} + P^{-1} r^{(k)}
\end{aligned}
$$

## Jacobi and JOR methods

For the Jacobi method $P = D$ and $N = D - A = E + F$ where $E$ and $F$ are the opposite of the lower and upper triangular parts of $A$, respectively. Component wise the iteration is as follows:

$$
x_i^{(k+1)} = \frac{1}{a_{ii}} [b_i - \sum_{j \neq i} a_{ij} x_j^{(k)}]
$$

The iteration matrix is

$$
B = D^{-1}(E + F) = I - D^{-1} A
$$

A relaxation parameter $\omega$ can be added to the Jacobi method to obtain the over-relaxation method (or JOR)

$$
x_i^{(k+1)} = \frac{\omega}{a_{ii}} [b_i - \sum_{j \neq i} a_{ij} x_j^{(k)}] + (1 - \omega)x_i^{(k)}
$$

and the iteration matrix is

$$
B_{J\omega} = \omega B_J + (1 - \omega)I
$$


## Gauss-Seidel and SOR methods

The Gauss-Seidel method uses triangular parts of the initial matrix for iterations such that $P = D - E$ and $N = F$. Component wise the iteration is expanded as

$$
x_i^{(k+1)} = \frac{1}{a_{ii}} [b_i - \sum_{j = 1}^{i-1} a_{ij} x_j^{(k+1)} - \sum_{j = 1}^{i-1} a_{ij} x_j^{(k)}]
$$

The iteration matrix is $B = (D - E)^{-1} F$

Same as for the Jacobi method a relaxation parameter $\omega$ can be added to yield the successive over-relaxation method (SOR)

$$
x_i^{(k+1)} = \frac{\omega}{a_{ii}} [b_i - \sum_{j = 1}^{i-1} a_{ij} x_j^{(k+1)} - \sum_{j = 1}^{i-1} a_{ij} x_j^{(k)}] + (1 - \omega) x_i^{(k)}
$$

for which the iteration matrix is

$$
B = (I - \omega D^{-1}E)^{-1}[(\omega D^{-1}F + (1 - \omega)I)x^{(k)} + b]
$$

## Convergence results for Jacobi, Gauss-Seidel, JOR and SOR

If $A$ is a strictly diagonally dominant matrix by rows, then the Jacobi and Gauss-Seidel methods are convergent.

If the Jacobi method is convergent then JOR method is convergent for $0 < \omega \leq 1$.
