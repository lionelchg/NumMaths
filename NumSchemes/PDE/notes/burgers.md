# Burgers equation notes

Let us consider the Burgers equation

$$
u_t + u u_x = 0
$$

So that the solution is

$$
u(x, t) = u_0(x_0)
$$

where

$$
x = x_0 + u_0(x_0) t
$$

## Smooth solution to shock

We consider a sin profile leading to shock at $t_s$ in the $[0, 1]$ domain

$$
u_0(x) = \frac{1}{2\pi t_s} \sin(2 \pi x)
$$

## Smooth rarefaction

An increasing profile leading to a rarefaction

$$
u_0(x) = A \tanh(k(x - 0.5))
$$

## Compression turning into a shock

Compression turning into a shock at $(x_s, t_s) = (0.75, 0.5)$

$$
\begin{aligned}
u_0(x) =
\begin{cases}
3/2 \quad x \leq 0 \\
3/2 - 2x \quad 0 < x < 1\\
-1/2 \quad x \geq 1
\end{cases}
\end{aligned}
$$

The solution before shock is

$$
\begin{aligned}
u_0(x) =
\begin{cases}
3/2 \quad x \leq 1.5 t \\
-1/2 \quad x \geq 1 - 0.5 t\\
\frac{3 - 4x}{2 - 4t}
\end{cases}
\end{aligned}
$$

The solution after the shock is

$$
\begin{aligned}
u_0(x) =
\begin{cases}
3/2 \quad x \leq x_s + 0.5 t \\
-1/2 \quad x \geq x_s + 0.5 t\\
\end{cases}
\end{aligned}
$$