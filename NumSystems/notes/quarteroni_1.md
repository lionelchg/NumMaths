# Chapter 1: Foundations of Matrix Analysis

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