# Numerical Computation of Internal and External Flows notes

## Chapter 8: General properties and High-Resolution Numerical Schemes

### Two level explicit schemes

We consider a general two level explicit scheme of the general form:

$$
    u_i^{n + 1} = \sum_j b_j u_{i+j}^n
$$

The range of $j$ is called the support of the schemes and is separated between $ju$ upwind points and $jd$ downwind points so that the total number of support points is $M=ju + jd + 1$.

We recall the first order upwind scheme and second order Lax-Wendroff schemes for the linear convection equation:

$$
\begin{aligned}
    &\mrm{FOU} \qquad u_i^{n+1} = u_i^n - \sig(u_i^n - u_{i-1}^n) \\
    &\mrm{LW} \qquad u_i^{n+1} = u_i^n - \frac{\sig}{2}(u_{i+1}^n - u_{i-1}^n) + \frac{\sig^2}{2}(u_{i+1}^n + u_{i-1}^n - 2u_i^n)
\end{aligned}
$$

where $(b_{-1}, b_0) = (\sig, 1 - \sig)$ for the FOU and $(b_{-1}, b_0, b_1) = (\sig(\sig + 1)/2, 1 - \sig^2, \sig(\sig - 1)/2)$ for LW.

We will now search for a general formulation of two level explicit schemes. Starting from Taylor expansions in time and space:

$$
\begin{aligned}
    u_i^{n+1} &= \sum_{m = 0}^{+\infty} \frac{\Delta t^m}{m!} \left(\pmdv{u}{t}{m}\right) \\
    u_{i+j}^n &= \sum_{m = 0}^{+\infty} \frac{(j \Delta x)^m}{m!} \left(\pmdv{u}{x}{m}\right) \\
\end{aligned}
$$

The first requirement is that a constant function should be solution of the problem:

$$
    \sum_j b_j = 1
$$

Reinjecting the Taylor expansions into the initial formula:

$$
    \Delta t \pdv{u}{t} + \sum_{m=2}^{+\infty} \frac{\Delta t^m}{m!} \left(\pmdv{u}{t}{m}\right) = \sum_j j b_j \Delta x \pdv{u}{x} + \sum_j b_j \sum_{m = 2}^{+\infty} \frac{(j \Delta x)^m}{m!} \left(\pmdv{u}{x}{m}\right)
$$

For the linear convection equation:

$$
    \pdv{u}{t} = - a \pdv{u}{x} \implies - a \Delta t \pdv{u}{x} + \ldots = \sum_j j b_j \Delta x \pdv{u}{x} + \ldots \implies \sum_j j \cdot b_j = - \sig 
$$

This yields:

$$
\begin{aligned}
    &\pdv{u}{t} = -a \pdv{u}{x} + \frac{1}{\Delta t}\sum_j b_j \sum_{m = 2}^{+\infty} \frac{(j \Delta x)^m}{m!} \left(\pmdv{u}{x}{m}\right) - \sum_{m=2}^{+\infty} \frac{\Delta t^{m-1}}{m!} \left(\pmdv{u}{t}{m}\right) \\
\implies &\pdv{}{t} = -a \pdv{}{x} + \frac{1}{\Delta t}\sum_j b_j\sum_{m = 2}^{+\infty} \frac{(j \Delta x)^m}{m!} \left(\pmdv{}{x}{m}\right) - \sum_{m=2}^{+\infty} \frac{\Delta t^{m-1}}{m!} \left(\pmdv{}{t}{m}\right)
\end{aligned}
$$

Let us recall a simple formula where $y$ and $z$ are first order terms compared to $x$:

$$
    (x + y + z)^m = x^m\left(1 + \frac{y}{x} + \frac{z}{x} \right)^m = x^m + mx^{m-1} y + mx^{m-1} z + \mrm{HOT}
$$

Developing the derivative operator using the formula above:

$$
    \pdv{u}{t} + a \pdv{u}{x} = \sum_{m = 2}^{+\infty} \alpha_m \Delta x^{m-1} \left(\pmdv{u}{x}{m}\right) - \sum_{m=2}^{+\infty} \frac{(-\sig)^{m-1}}{(m-1)!} \sum_{n=2}^{+\infty} \alpha_n \Delta x^{m+n-2} \left(\pmdv{u}{x}{n+m-1}\right)
$$

where 

$$
    \alpha_p = \frac{1}{p!} \left[\sum_j b_j j^p - (-\sig)^p\right] \frac{\Delta x}{\Delta t}
$$

Then requiring the scheme to be of order $p$ yields the condition:

$$
    \pdv{u}{t} + a \pdv{u}{x} = \mc{O}(\Delta x^p) \implies \sum_j b_j j^m = (-\sig)^m \:\: \mrm{for}\: m \in \intint{0}{p}
$$

When the scheme is of order $p$ then the numerical scheme writes:

$$
    \pdv{u}{t} + a \pdv{u}{x} = \sum_{m=p}^{+\infty} a_{m+1} \Delta x^m \left(\pmdv{u}{x}{m+1}\right)
$$

and the $a_j$ are linked with the $\alpha_j$:

$$
\begin{aligned}
    a_{p+1} &= \alpha_{p+1} \\
    a_{p+2} &= \alpha_{p+2} + \sig \alpha_{p+1} \\
    a_{p+3} &= \alpha_{p+3} + \sig \alpha_{p+2} - \frac{\sig^2}{2}\alpha_{p+1} \\
    a_{p+4} &= \alpha_{p+4} + \sig \alpha_{p+3} - \frac{\sig^2}{2}\alpha_{p+2} + \frac{\sig^3}{6}\alpha_{p+1}\\
\end{aligned}
$$

### One parameter family schemes on the support $(i - 2, i - 1, i, i + 1)$

We consider second order schemes on the support $(i - 2, i - 1, i, i + 1)$. Thus the $b_j$ verify the relations:

$$
\begin{aligned}
b_{-2} + b_{-1} + b_0 + b_1 &= 1 \\
-2b_{-2} - b_{-1} + b_1 &= -\sig \\
4b_{-2} + b_{-1} + b_1 &= \sig^2
\end{aligned}
$$

Solving the system:

$$
\begin{aligned}
b_{-2} &= - \gamma\\
b_{-1} &= \frac{\sig(\sig + 1)}{2} + 3 \gamma\\
b_0 &= 1 - \sig^2 - 3 \gamma\\
b_1 &= \frac{\sig(\sig - 1)}{2} + \gamma
\end{aligned}
$$

Hence the scheme can be written as:

$$
    S = S_\mrm{LW} + \gamma H(-1, 3, -3, 1)
$$

Multiple schemes can then be derived:
$$
\begin{aligned}
    &S_\mrm{LW}(0, \frac{\sig(\sig + 1)}{2}, 1 - \sig^2, \frac{\sig(\sig - 1)}{2}) \\
    &S_\mrm{SOU}(\frac{\sig(\sig - 1)}{2}, \sig(2 - \sig), \frac{(1 - \sig)(2 - \sig)}{2}, 0) \\
    &S_\mrm{Fromm}(\frac{\sig(\sig - 1)}{4}, \frac{\sig(5 - \sig)}{4}, \frac{(1- \sig)(4 + \sig)}{4}, \frac{\sig(\sig - 1)}{4}) \\
    &S_\mrm{3}(\frac{\sig(\sig^2 - 1)}{6}, \frac{\sig(2 - \sig)(\sig + 1)}{2}, \frac{\sig(2 - \sig)(1 - \sig^2)}{2}, \frac{\sig(\sig - 1)(2 - \sig^2)}{6})
\end{aligned}
$$