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