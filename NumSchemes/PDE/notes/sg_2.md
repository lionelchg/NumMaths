# Multi-dimensional SG scheme


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

Applying the SG scheme in a multi-dimensional setting

$$
\begin{aligned}
\vb{F} &= \vb{V} u - D \nabla u \\
\vb{F}_{ij} \cdot \hat{\vb{ij}} &= \vb{V}_{ij} \cdot \hat{\vb{ij}} u - D \nabla u \cdot \hat{\vb{ij}}\\
\end{aligned}
$$

The SG scheme only gives tangentiel flux projection:

$$
S_\mrm{SG}(\vb{F}) = \vb{V}_{ij} \cdot \vb{S}_{ij}^{\parallel} \frac{u_j - e^\alpha u_i}{1 - e^\alpha} 
$$

where 

$$
\alpha = \frac{h_{ij} \vb{V}_{ij} \cdot \hat{\vb{ij}}}{D_{ij}}
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

On meshes where the perpendicular component is non-zero another scheme must be applied on top of the SG scheme. A solution seems to simply define the centered scheme for it:

$$
\vb{F}_{ij} \cdot \vb{S}_{ij}^\perp = \vb{V}_{ij} \frac{u_i + u_j}{2} \cdot \vb{S}_{ij}^\perp - D_{ij} \frac{\nabla u_i + \nabla u_j}{2} \cdot \vb{S}_{ij}^\perp
$$

The instabilities are mostly arising in zones where $\alpha \ll 1$ where the diffusion dominates and in that case the scheme reduces to it. Now is there an issue to have an upwind in one direction and central difference in another?