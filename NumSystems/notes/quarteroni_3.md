# Chapter 3 - Direct methods for the solution of linear systems

In this chapter we are interested in linear systems, we want to find $x$ such that

$$
Ax = b
$$

We restrict the study to square systems where $A$ is invertible so that the system admits a unique solution. Exact solution of such systems are given by Cramer's rule:

$$
x_j = \frac{\Delta_j}{\det(A)}
$$

where $\Delta_j$ is the determinant of $A$ where the $j$-th column is replaced by $b$. This formula can not be used in practice due to its $\mc{O}((n+1)!)$ complexity.

## Stability analysis

The condition number of a matrix is defined as follows

$$
K(A) = ||A|| \: ||A^{-1}||
$$

The condition number for the 2-norm is called the spectral condition number because

$$
K_2(A) = \frac{\sigma_1(A)}{\sigma_n(A)}
$$

The distance of matrix $A$ to the set of singular matrices is

$$
\mrm{dist}(A) = \min\{||\delta A||/||A||: \det(A + \delta A) = 0\} = 1 / K(A)
$$

Hence for invertible matrices:

$$
||\delta A|| \: ||A^{-1}|| < 1
$$

### Forward analysis

Given $\delta A$ and $\delta b$ we adjust $x$ such that

$$
(A + \delta A)(x + \delta x) = b + \delta b
$$

The following relation holds

$$
\frac{||\delta x||}{||x||} \leq \frac{K(A)}{1 - K(A)||\delta A||/||A||}\left(\frac{||\delta b||}{||b||} + \frac{||\delta A||}{||A||}\right)
$$

### Backward analysis

The numerical method yields a solution $\hat{x} = Cb$ where the matrix $C$, due to rounding errors, is an approximation of $A^{-1}$. Let $R = AC - I$. If $||R|| < 1$ then $A$ and $C$ are nonsingular and

$$
||A^{-1}|| \leq \frac{||C||}{1 - ||R||}, \frac{||R||}{||A||} \leq ||C - A^{-1}|| \leq \frac{||C||\:||R||}{1 - ||R||}
$$

### A posteriori analysis

Let us denote by $y$ an approximate solution. The error $e = y - x$ can be related to $y$, $C$ and the residual vector $r = b - Ay$

$$
||e|| \leq \frac{||r|| \: ||C||}{1 - ||R||}
$$

