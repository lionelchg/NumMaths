# Scharfetter Gummel scheme

## Derivation

We consider the general advection diffusion equation

$$
\pdv{u}{t} + \pdv{F}{x} = 0
$$

where

$$
F = V u - D \pdv{u}{x}
$$

Placing ourselves in a mesh between $x_i$ and $x_{i+1}$ the first order differential equation 

$$
\pdv{u}{x} - \frac{V_{i+1/2}}{D_{i+1/2}} u = - \frac{F_{i+1/2}}{D_{i+1/2}}
$$

is solved. The solution is

$$
u(x) = \left[u_i - \frac{h_iF}{D} \int_0^\xi \exp(-\alpha \xi') \mrm{d}\xi' \right] \exp(\alpha \xi)
$$

where

$$
\begin{aligned}
\xi &= \frac{x - x_i}{h_i} \\
\alpha &= \frac{h_i V_{i+1/2}}{D_{i+1/2}} \\
h_i &= x_{i+1} - x_i
\end{aligned}
$$

Evaluating the solution at $x_{i+1}$ gives

$$
F_{i+1/2} = \frac{V_{i+1/2}}{1 - e^\alpha}[u_{i+1} - e^\alpha u_i]
$$

## Asymptotic regimes

### $\alpha \gg 1$ 

The flux reduces to a first order upwind for the convection part

$$
F_{i+1/2} = V_{i+1/2}u_i
$$

### $\alpha \ll 1$

The flux reduces to centered difference for the diffusive flux and convective flux

$$
F_{i+1/2} = - \frac{D_{i+1/2}}{h_i} (u_{i+1} - u_i) + V_{i+1/2} \frac{u_{i+1} + u_i}{2} - \frac{\alpha}{12} V_{i+1/2}(u_{i+1} - u_i)
$$

## Linear velocity

Assuming a linear velocity profile between $x_i$ and $x_{i+1}$

$$
V(x) = V_i(1 + \frac{\Delta V_i}{V_i}\xi) = V_i(1 + 2\beta\xi)
$$

with $\beta = \Delta V_i / 2 V_i$

Solving the first order differential equation as above with the linear profile yields

$$
u(x) = \left[u_i - \frac{h_iF_{i+1/2}}{D_{i+1/2}} \int_0^\xi \exp(-f(\xi')) \mrm{d}\xi' \right] \exp(f(\xi))
$$

where

$$
\begin{aligned}
\xi &= \frac{x - x_i}{h_i} \\
f(\xi) &= \alpha(1 + \beta \xi) \xi \\
\alpha &= \frac{V_i h_i}{D_{i+1/2}}
\end{aligned}
$$

Assuming $|\alpha \beta| \ll 1$

$$
F_{i+1/2} = \frac{D_{i+1/2}}{h_i I_1}\left[u_i - \frac{\exp(-\alpha \xi)}{1 + \alpha \beta} u_{i+1}\right]
$$

with

$$
I_1 = \int_0^1 \exp(-\alpha \xi')(1 - \alpha \beta \xi'^2)\mrm{d}\xi'.
$$

This integral can be calculated exactly

$$
I_1 = \frac{\alpha^2 - 2 \alpha \beta}{\alpha^3} + \frac{e^{-\alpha}}{\alpha^3}(\alpha^3 \beta - \alpha^2 + 2 \alpha \beta + 2 \alpha^2 \beta)
$$

## Virtual nodes

To satisfy $|\alpha\beta| \ll 1$, virtual nodes are defined with a parameter $\veps$

$$
h_v = \sqrt{\frac{\veps 2 D_{i+1/2}h_i}{|\Delta V_i|}} 
$$

The coordinates and values of velocities are defined as

$$
\begin{aligned}
x_{L, R} &= x_{i+1/2} \mp \frac{h_v}{2} \\
V_{L, R} &= V_{i+1/2} \mp \frac{\Delta V_v}{2} \\
\Delta V_v &= \frac{h_v}{h_i} \Delta V_i
\end{aligned}
$$

Define the interpolated scalar values $u_{L, R}$ with exponential interpolation

$$
u(x) = u_i \exp(a(x - x_i)),\quad a = \frac{1}{h_i}\ln(\frac{n_{i+1}}{n_i})
$$