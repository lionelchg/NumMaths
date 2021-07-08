# Numerical Computation of Internal and External Flows notes

## Chapter 8: General properties and High-Resolution Numerical Schemes

### Monotonicity of Numerical Schemes

#### Monotonicity conditions

Taking a two-step numerical scheme:

$$
    u_i^{n+1} = \sum_j b_j u_{i+j}^n
$$

A scheme is said to be monotone if the new solution value at time $(n+1)$, $u_i^{n+1}$, does not reach values outside the range covered by $u_{i+j}^n$. 

Then

$$
u_\mrm{min}^n \leq u_{i+j}^n \leq u_\mrm{max}^n \implies u_\mrm{min}^n \leq u_{i}^{n+1} \leq u_\mrm{max}^n
$$

only if $b_j \geq 0$ for all $j$ in the support of the scheme.

An equivalent way of looking at monotonicity is to see it as the requirement that a local minimum cannot decrease and a local maximum cannot increase when going from one timestep to the other. Rewriting the previous two-step numerical scheme with $\sum_j b_j = 1$:

$$
    u_i^{n+1} = u_i^n + \sum_j b_j (u_{i+j}^n - u_i^n)
$$

If $b_j \geq 0$ is satisfied then $u_i^{n+1} \geq u_i^n$ for a local minimum and $u_i^{n+1} \leq u_i^n$ for a local maximum.

#### Semi-discretized schemes or method of lines

To generalize the results of the previous section, we consider semi-discretized schemes of the form:

$$
\dv{u_i}{t} = \sum_j \beta_j u_{i+j}
$$

For a constant value to be solution of the scheme:

$$
\sum_j \beta_j = 0
$$

Then:

$$
\dv{u_i}{t} = \sum_{j\neq 0} \beta_j (u_{i+j} - u_i)
$$

Hence the monotonicity, positivity or Local Extremum Diminishing (LED) condition is:

$$
\forall j \neq 0 \: \beta_j \geq 0
$$

The dominating problem for monotonicity of a scheme comes from the convection term and focus should be made on this term.

#### Godunov's theorem

***Theorem***: All linear monotone schemes for the convection equation are necessarily first order accurate.

The proof is only made for the two-step explicit schemes. The requirement for a second order scheme is:

$$
\begin{aligned}
\sum_j j b_j &= - \sig \\
\sum_j j^2 b_j &= \sig^2 \\
\end{aligned}
$$

Using Schwart's inequality for a first order monotone scheme:

$$
\begin{aligned}
\sig^2 = \left(\sum_j j b_j\right)^2 = \left(\sum_j j \sqrt{b_j} \sqrt{b_j}\right)^2 \leq \sum_j j^2 b_j
\end{aligned}
$$

Moreover since the two vectors $(j\sqrt{b_j})$ and $(\sqrt{b_j})$ are non colinear (except for $b_j = 0$ which is not possible) then the inequality is strict and $\sig^2 < \sum_j j^2 b_j$ finishing the proof.

#### High-resolution schemes and the concept of limiters

1. Select a first-order monotone scheme (usually the upwind scheme), as reference. Express the high order scheme as the monotone scheme plus additional terms
2. Multiply the additional terms by a limiting funciton $\Psi(r_i)$, expressed as a function of ratios of successive gradients.
3. Express the monotonicity conditions to derive conditions on the limiters

For example for a second order upwind scheme:
$$
\begin{aligned}
\dv{u_i}{t} &= - \frac{a}{2 \Delta x}(3u_i - 4 u_{i-1} + u_{i-2}) \\
&= - \frac{a}{2\Delta x}(-4(u_{i-1} - u_i) + (u_{i-2} - u_i)) \\
&= - \frac{a}{\Delta x}((u_i - u_{i-1}) - \frac{a}{2\Delta x}(-2(u_{i-1} - u_i) + (u_{i-2} - u_i)) \\
&= - \frac{a}{\Delta x}((u_i - u_{i-1}) - \frac{a}{2\Delta x}(-(u_{i-1} - u_i) + (u_{i-2} - u_{i-1})) \\
\end{aligned}
$$

We introduce the ratio of successive gradient:
$$
r_i = \frac{u_{i+1} - u_i}{u_{i} - u_{i-1}}
$$

Multiplying the additional terms:

$$
\begin{aligned}
\dv{u_i}{t} &= - \frac{a}{\Delta x}((u_i - u_{i-1}) - \frac{a}{\Delta x}(\frac{1}{2}(u_i - u_{i-1}) - \frac{1}{2}(u_{i-1} - u_{i-2})) \\
&= - \frac{a}{\Delta x}((u_i - u_{i-1}) - \frac{a}{\Delta x}(\frac{1}{2}\Psi(r_i)(u_i - u_{i-1}) - \frac{1}{2}\Psi(r_{i-1})(u_{i-1} - u_{i-2})) \\
&= - \frac{a}{\Delta x}\left[1 - \frac{1}{2}\Psi(r_i) + \frac{1}{2}\Psi(r_{i-1})/r_{i-1}\right](u_i - u_{i-1})
\end{aligned}
$$

Requiring the terms in the brackets to be positive yields:

$$
    0 \leq \Psi(r) \leq \min(2r, 2)
$$

For the accuracy of a linear solution:

$$
    \Psi(1) = 1
$$

For the SOU scheme as base:

$$
u_i^{n+1} = u_i^n - \sigma\left[1 + \frac{1}{2}(1-\sigma)(\Psi(r_i) - \Psi(r_{i-1})/r_{i-1})\right]
$$

We recover the SOU scheme for $\Psi = 1$ and the LW scheme for $\Psi = r$.

Various limiters can then be defined:

$$
\begin{aligned}
&\mrm{Van Leer} \qquad \Psi(r) = \frac{r + |r|}{1 + r} \\
&\mrm{Min-mod} \qquad \Psi(r) = \mrm{min-mod}(r, 1) \\
&\mrm{Superbee} \qquad \Psi(r) = \max(0, \min(2r, 1), \min(r, 2)) \\
&\beta-\mrm{limiters} \qquad \Psi(r) = \max(0, \min(\beta r, 1), \min(r, \beta)) \\
&\mrm{Osher limiter} \qquad \Psi(r) = \max(0, \min(r, 2)) \\
&\alpha-\mrm{limiters} \qquad \Psi(r) = \max(0, \min(2r, \alpha r + 1 - \alpha, 2))
\end{aligned}
$$