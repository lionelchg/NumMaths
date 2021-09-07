# Chapter 1: Foundations of Matrix Analysis

## General definitions

The maximum module of the eigenvalues of a matrix $A$ is called the spectral radius $\rho(A)$:

$$
\rho(A) = \max_{\lambda \in \mrm{Sp}(A)} |\lambda|
$$

## Similarity transformations

### Schur decomposition

Given $A \in \mc{M}_n(\mathbb{C})$ there exists $U$ unitary such that

$$
U^{-1} A U = U^H A U =
\begin{bmatrix}
\lambda_1 & \cdots & & b_{1n} \\
0 & \lambda_2 & & b_{2n} \\
\vdots & & \ddots & \vdots \\
0 & \cdots & 0 & \lambda_n
\end{bmatrix}
$$

where $\lambda_i$ are the eigenvalues of the matrix.

### Canonical Jordan form

Let $A$ be any square matrix. There exists a nonsingular matrix $X$ that transforms $A$ into a block diagonal matrix $J$

$$
X^{-1} A X = J = \mrm{diag}(J_{k_1}(\lambda_1), \ldots, J_{k_l}(\lambda_l))
$$

which is called canonical Jordan form where $\lambda_j$ are the eigenvalues of $A$ and $J_k(\lambda) \in \mc{M}_k(\mathbb{C})$ with $J_k(\lambda) = 1$ if $k = 1$ and 

$$
J_k(\lambda) = 
\begin{bmatrix}
\lambda & 1 & 0 & \ldots & 0 \\
0 & \lambda & 1 & \ldots & 0 \\
0 & 0 & \lambda &  & 0 \\
\vdots & \vdots & & \ddots & 1 \\
0 & 0 &  & \ldots & \lambda \\
\end{bmatrix}
$$

otherwise.

### Singular Value Decomposition

Let $A \in \mc{M}_{mn}(\CC)$. There exists two unitary matrices $U \in \mc{M}_m(\CC)$ and $V \in \mc{M}_n(\CC)$ such that

$$
U^H A V = \Sigma = \mrm{diag}(\sigma_1, \ldots, \sigma_p) \in \mc{M}_{mn}(\CC)
$$
$\Sigma$ is called the Singular Value decomposition of $A$ and $\sigma_1 > \sigma_2 > \ldots > \sigma_p$ are called the singular values of $A$.

## Matrix Norms

### Definitions

A matrix norm $||.||$ is said to be compatible or consistent with a vector norm $||.||$ if:

$$
\begin{equation}
\forall x \in \R^n \ ||Ax|| \leq ||A|| \: ||x||
\end{equation}
$$

A matrix norm is said to be submultiplicative if

$$
||AB|| \leq ||A|| \, ||B||
$$

The Frobenius norm

$$
||A||_F = \mathrm{tr}(AA^H)
$$
is compatible with the 2-norm.

Let $||.||$ be a vector norm. The induced or natural matrix norm is:

$$
||A|| = \sup_{x \neq 0} \frac{||Ax||}{||x||}
$$

The 1-norm and infinity-norm as easily computable since:

$$
\begin{align}
||A||_1 &= \max_j \sum_i |a_{ij}| \\
||A||_\infty &= \max_i \sum_j |a_{ij}|
\end{align}
$$

The 2-norm computation is much more expensive. Denoting by $\sigma_1(A)$ the largest singular value of $A$:

$$
||A||_2 = \sqrt{\rho(A^HA)} = \sqrt{\rho(AA^H)} = \sigma_1(A)
$$

In particular if $A$ is hermitian:

$$
||A||_2 = \rho(A)
$$

### Properties

Let $||.||$ be a consistent matrix norm, then:

$$
\forall A \in \mathbb{C}^{n \times n} \ \rho(A) \leq ||A||
$$

More precisely:

$$
\rho(A) = \inf_{||.||} ||A||
$$

A sequence of matrix $A_k$ is said to converge to a matrix $A$ if

$$
\lim_{k\to +\infty} ||A_k - A|| = 0
$$

The following equivalences hold

$$
\lim_k A^k = 0 \iff \rho(A) < 1 \iff \sum_k A^k = (I - A)^{-1}
$$

Moreover if $\rho(A) < 1$:

$$
\frac{1}{1 + ||A||} \leq ||(I - A)^{-1}|| \leq \frac{1}{1 - ||A||}
$$