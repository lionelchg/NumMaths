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
x_i^{(k+1)} = \frac{1}{a_{ii}} [b_i - \sum_{j = 1}^{i-1} a_{ij} x_j^{(k+1)} - \sum_{j = i+1}^{n} a_{ij} x_j^{(k)}]
$$

The iteration matrix is $B = (D - E)^{-1} F$

Same as for the Jacobi method a relaxation parameter $\omega$ can be added to yield the successive over-relaxation method (SOR)

$$
x_i^{(k+1)} = \frac{\omega}{a_{ii}} [b_i - \sum_{j = 1}^{i-1} a_{ij} x_j^{(k+1)} - \sum_{j = i+1}^{n} a_{ij} x_j^{(k)}] + (1 - \omega) x_i^{(k)}
$$

for which the iteration matrix is

$$
B = (I - \omega D^{-1}E)^{-1}[(\omega D^{-1}F + (1 - \omega)I)x^{(k)} + b]
$$

## Convergence results for Jacobi, Gauss-Seidel, JOR and SOR

If $A$ is a strictly diagonally dominant matrix by rows, then the Jacobi and Gauss-Seidel methods are convergent.

If the Jacobi method is convergent then JOR method is convergent for $0 < \omega \leq 1$.

## Stationary and Non-stationary methods

The non-stationary Richardson method is the following

$$
x^{(k+1)} = x^{(k)} + \alpha_k P^{-1} r^{(k)}
$$

$\alpha_k = \alpha$ corresponds to the stationary Richardson method. The associated iteration matrix is

$$
B = I - \alpha_k P^{-1}A
$$

For computations the following steps are performed

1. Solve $Pz^{(k)} = r^{(k)}$
2. Compute $\alpha_k$
3. Update solution $x^{(k+1)} = x^{(k)} + \alpha_k z^{(k)}$
4. Update residual $r^{(k+1)} = r^{(k)} - \alpha_k A z^{(k)}$

For any non-singular matrix $P$ the stationary Richardson method is convergent if

$$
\frac{2 \mrm{Re}\lambda_i}{\alpha |\lambda_i|^2} > 1
$$

where $\lambda_i$ are the eigenvalues of $P^{-1}A$.

## Preconditioning matrices

To accelerate the convergence procedure, preconditioning can be applied in different forms:

1. Left preconditioning $P^{-1} A x = P^{-1} b$
2. Right preconditioning $AP^{-1} y = b$ and $y = Px$
3. Centered preconditioning $P_L^{-1}AP_R^{-1} y = b$ and $y = P_Rx$

## The Gradient method

For a positive definite matrix $A$ we define the energy of the linear system

$$
E(y) = \frac{1}{2} y^T A y - y^Tb
$$

The gradient is $\nabla E = Ay - b$ so that for the solution of the system $x$ we have $\nabla E = 0$ and from Taylor expansion in multiple dimensions

$$
E(y) = E(x + (y-x)) = E(x) + \frac{1}{2} ||y - x||_A^2
$$

We use a Richardson method and choose the direction of update as the direction of steepest gradient descent

$$
x^{(k+1)} = x^{(k)} + \alpha_k r^{(k)}
$$

The optimal parameter is

$$
\alpha_k = \frac{r^{(k)T}r^{(k)}}{r^{(k)T}Ar^{(k)}}
$$
