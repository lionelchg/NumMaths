# Chapter 2 - Principles of Numerical Mathematics

## Well-posedness and condition number of a problem

Let us consider the problem to find $x$ such that

$$
F(x, d) = 0
$$

where $d$ is the set of data which the solution depends on and $F$ is the functional relationship between $x$ and $d$. We restrict the study to direct problems where $F$ and $d$ are given and $x$ is the unkwown.

The problem is well-posed if it admits a unique solution which depends continuously on the data. Continuous dependence on the data means that small perturbations of $d$ leads to small perturbations in the solution $x$. For a given $\delta d$, the subsequent change $\delta x$ such that

$$
F(x + \delta x, d + \delta d) = 0
$$

satisfies

$$
\forall \eta > 0 \ \exists K(\eta, d) \ ||\delta d|| < \eta \implies ||\delta x|| < K(\eta, d) ||\delta d||
$$

Following that, the relative $K(d)$ and absolute $K_\mrm{abs}(d)$ condition numbers are defined

$$
\begin{align}
K(d) &= \sup_{\delta d \in D} \frac{||\delta x||/||x||}{||\delta d||/||d||}\\
K_\mrm{abs}(d) &= \sup_{\delta d \in D} \frac{||\delta x||}{||\delta d||}
\end{align}
$$

where $D$ is a neighborhood of the origin for which the perturbed system still makes sense. Introducing the resolvent $G$ for a well-posed problem

$$
x = G(d)
$$

The condition number are simplified

$$
\begin{align}
K(d) &= \frac{||G'(d)|| \, ||d||}{||G(d)||}\\
K_\mrm{abs}(d) &= ||G'(d)||
\end{align}
$$

## Stability and convergence of numerical methods

### Stability

We suppose the problem $F(x, d) = 0$ to be well-posed. A numerical method for solving the problem will consist in a sequence of approximate problems

$$
F_n(x_n, d_n) = 0 \qquad n \geq 1
$$

The expectation is that $x_n \to x$, $d_n \to x$ and $F_n$ approximates $F$ as $n \to +\infty$.

The numerical method is said to be consistent if

$$
F_n(x, d) = F_n(x, d) - F(x, d) \to 0 \qquad \text{when } n \to +\infty
$$

The numerical method is said to be strongly consistent if $F_n(x, d) = 0$ for any value of $n$.

We introduce the same definition of well-posedness of a numerical method: for any fixed value of $n$, there exists a unique solution $x_n$ given data $d_n$ that depends continuously on the data

$$
\forall \eta > 0 \ \exists K_n(\eta, d) \ ||\delta d_n|| < \eta \implies ||\delta x_n|| < K_n(\eta, d_n) ||\delta d_n||
$$

Following the same definitions as in the first part

$$
\begin{align}
K_n(d_n) &= \sup_{\delta d_n \in D_n} \frac{||\delta x_n||/||x_n||}{||\delta d_n||/||d_n||}\\
K_{\mrm{abs}, n}(d) &= \sup_{\delta d_n \in D_n} \frac{||\delta x_n||}{||\delta d_n||}
\end{align}
$$

and then

$$
\begin{align}
K^\mrm{num}(d_n) &= \lim_{k \to +\infty} \sup_{n \geq k} K_n(d_n) \\
K^\mrm{num}_\mrm{abs}(d_n) &= \lim_{k \to +\infty} \sup_{n \geq k} K_{\mrm{abs}, n}(d_n)
\end{align}
$$

We call $K^\mrm{num}(d_n)$ the relative asymptotic condition number of the numerical method and $K^\mrm{num}_\mrm{abs}(d_n)$ the asbolute asymptotic condition number of the numerical method, corresponding to the data $d_n$.

Introducing the numerical resolvent $G_n$ for a well-posed problem

$$
x_n = G_n(d_n)
$$

The condition number are simplified

$$
\begin{align}
K_n(d_n) &= \frac{||G_n'(d_n)|| \, ||d_n||}{||G_n(d_n)||}\\
K_{\mrm{abs}, n}(d) &= ||G_n'(d_n)||
\end{align}
$$

### Convergence

The numerical method is said to be convergent if

$$
\begin{aligned}
&\forall \veps > 0 \ \exists n_0(\veps) \ \exists \delta(n_0, \veps) > 0: \\
&\forall n > n_0(\veps) \ \forall ||\delta d_n|| < \delta(n_0, \veps) \implies ||x(d) - x_n(d + \delta d_n)|| < \veps
\end{aligned}
$$
where $d$ is an admissible data for the problem, $x(d)$ is the corresponding solution and $x_n(d_n + \delta d_n)$ is the solution of the numerical problem with data $d_n + \delta d_n$.

The link between the two is the Lax-Richtmeyer theorem: for a consistent numerical method, stability is equivalent to convergence.

## Stability analysis

### A priori analysis

1. *Forward analysis*: it provides a bound to the variations $||\delta x_n||$ on the solution due to both perturbations in the data and to errors that are intrinsic to the numerical method
2. *Backward analysis*: given a certain computed solution $\hat{x}_n$, backward analysis looks for the perturbation on the data $\delta d_n$ such that $F_n(\hat{x}_n, d_n + \delta d_n) = 0$.

### A posteriori analysis

The a posteriori error aims at evaluating the error $\hat{x}_n - x$ as a function of the residual $r_n = F(\hat{x}_n, d)$ by means of constants that are called stability factors.
