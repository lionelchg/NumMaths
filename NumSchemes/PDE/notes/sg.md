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

## General scheme

For a general advection diffusion problem

$$
\pdv{u}{t} + \nabla \cdot \vb{F} = 0
$$

where

$$
\vb{F} = \vb{V} u - D \nabla u = \vb{F}_C + \vb{F}_D
$$

Integrating on the nodal volume in a vertex centered formulation

$$
\pdv{u_i}{t} + \frac{1}{V_i}\left[\int_{\partial \dot{V}_i} \vb{F}\cdot\vb{n}\mrm{d}S + \int_{\partial V_i^\Gamma} \vb{F}\cdot\vb{n}\mrm{d}S\right] = 0
$$

The scheme applies to the interior nodes so that

$$
\int_{\partial \dot{V}_i} \vb{F}\cdot\vb{n}\mrm{d}S = \sum_{e \in E(i)} \sum_{f \in K_e \cap \partial \dot{V}_i} \int_f \vb{F} \cdot \vb{n} \mrm{d}S
$$

Let us denote by $ij$ the selected edges so that the scheme is:

$$
\int_f \vb{F} \cdot \vb{n} \mrm{d}S = f(\vb{F}_i, \vb{F}_j)
$$

### Common schemes

$$
\begin{aligned}
S_\mrm{FOU}(\vb{F}_C) &= \max(0, \vb{S}_{ij} \cdot \vb{V}_{ij}) u_i + \min(0, \vb{S}_{ij} \cdot \vb{V}_{ij}) u_j \\
S_\mrm{CD}(\vb{F}_D) &= - D_{ij} \frac{\nabla u_i + \nabla u_j}{2} \cdot \vb{S}_{ij} \\
S_\mrm{CD}(\vb{F}_C) &= \vb{V}_{ij} \frac{u_i + u_j}{2} \cdot \vb{S}_{ij}
\end{aligned}
$$

The SG scheme only gives tangentiel flux projection:

$$
S_\mrm{SG}(\vb{F}) = \vb{V}_{ij} \cdot \vb{S}_{ij}^{\parallel} \frac{u_j - e^\alpha u_i}{1 - e^\alpha} 
$$

where 

$$
\alpha = \frac{\vb{V}_{ij} \cdot \vb{S}_{ij}^\parallel h_{ij}}{D_{ij}}
$$

### Flux splitting

Total flux can be splitted following:

$$
\int_f \vb{F} \cdot \vb{n} \mrm{d}S = \vb{F}_{ij} \cdot \vb{S}_{ij}
$$

Defining the tangential vector

$$
\hat{\vb{ij}} = \frac{\vb{ij}}{ij}
$$

The surface normal is splitted:

$$
\vb{S}_{ij} = \vb{S}_{ij}^\parallel + \vb{S}_{ij}^\perp
$$

So that total flux can be written

$$
\int_f \vb{F} \cdot \vb{n} \mrm{d}S = \vb{F}_{ij} \cdot \vb{S}_{ij}^\parallel + \vb{F}_{ij} \cdot \vb{S}_{ij}^\perp
$$