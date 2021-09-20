# Numerical Computation of Internal and External Flows notes

## Chapter 7: Consistency, Stability and Error Analysis of Numerical Schemes

### Basic concepts and definitions

#### Consistency

A numerical scheme $N(U_i^n) = 0$ is said to be consistent if it tends to the differential equation $D(U) = 0$ when time and space steps tend to zero:

$$
\lim_{\Delta x \to 0, \Delta t \to 0} N(U_i^n) = D(U_i^n)
$$

#### Stability

Stability is the requirement that all errors must remain bounded when the iteration process advances. Denoting by $u_i^n$ the computed solution and $\bar{u}_i^n$ the exact solution of the numerical scheme, we define the error as $\bar{\veps}_i^n = u_i^n - \bar{u}_i^n$ and require:

$$
    \lim_{n\to +\infty} |\bar{\veps}_i^n| \leq K 
$$
at fixed $\Delta t$ with $K$ independent of $n$.

#### Convergence

Convergence is a condition on the numerical solution. When time and space steps tend to zero, the numerical solution $u_i^n$ must converge to the exact solution $\tilde{u}_i^n$ of the differential equation. Defining the error as:

$$
    \tilde{\veps}_i^n = u_i^n - \tilde{u}_i^n
$$

we require that:

$$
    \lim_{\Delta x \to 0, \Delta t \to 0}\max_{n \in \intint{0}{T/\Delta t}, i \in \intint{0}{L_x/\Delta x}} |\tilde{\veps}_i^n| = 0
$$

## Truncation error

The truncation is obtained from

$$
N(u_i^n) - D(u_i^n) = \veps_T
$$

So that the equivalent differential equation is

$$
D(\bar{u}_i^n) = - \bar{\veps}_T
$$

### Von Neumann analysis

#### Fourier decomposition of the solution

Let us consider a domain of length $L$ in the $x$-axis with $N + 1$ points and spacing $\Delta x$. Reflecting this domain with respect to the origin the minimum and maximum wavelength and wave numbers are:

$$
\begin{aligned}
    &\lambda_\mrm{max} = 2L \implies k_\mrm{min} = \frac{\pi}{L} \\
    &\lambda_\mrm{min} = 2\Delta x \implies k_\mrm{max} = \frac{\pi}{\Delta x} \\
\end{aligned}
$$

Hence
$$
    k_j = \frac{j\pi}{L}, \: j \in \intint{-N}{N}
$$

And we can decompose the solution into spatial modes:

$$
    u_i^n = \sum_{j=-N}^N V_j^n e^{\I k_j x_i} = \sum_{j=-N}^N V_j^n e^{\I i \cdot k_j \Delta x} = \sum_{j=-N}^N V_j^n e^{\I i \cdot \phi_j}
$$

where $\phi_j = k_j \Delta x \in [-\pi, \pi]$ is the phase angle which remaps spatial frequencies into the range $[-\pi, \pi]$.

#### Amplification factor

The Von Neumann stability condition requires that the amplitudes do not grow in time, we define the amplification factor:

$$
    G = \frac{V^{n+1}}{V^n}
$$

and require that

$$
    |G(\phi_j)| \leq 1 \qquad \forall j \in \intint{-N}{N} 
$$

In practice for a numerical scheme we set:

$$
    u_{i+m}^{n+k} = V^{n+k} e^{\I(i+m)\phi}
$$

and find the amplification factor depending on $\phi$.

### Spectral analysis of numerical schemes

For the exact solution of the differential equation $\tu$ let us consider its mode of wave number $k$ where $\tilde{\omega} = \tilde{\omega}(k)$ is the exact dispersion relation:

$$
    \tu_k(x, t) = \hat{V}(k) e^{\I(kx - \tilde{\omega}t)}
$$

where

$$
    \hat{V}(k) = \frac{1}{2L} \int_{-L}^{L} u(x, 0)e^{-\I kx} dx
$$

Considering now the computed solution with an approximate dispersion relation $\omega = \omega(k)$:

$$
    (u_i^n)_k = \hat{V}(k) e^{\I kx_i} e^{-\I \omega t^n}
$$

We can then identify the exact and numerical amplification factors:

$$
\begin{aligned}
    \tilde{G} &= e^{-\I \tilde{\omega} \Delta  t} = |\tilde{{G}}| e^{-\I\tilde{\Phi}} \\
    G &= e^{-\I \omega \Delta  t} = |G| e^{-\I\Phi} 
\end{aligned}
$$

From these two amplification factors the diffusion and dispersion errors are:

$$
\begin{aligned}
    \veps_D &= \frac{|G|}{|\tilde{G}|} \\
    \veps_\phi &= \frac{\Phi}{\tilde{\Phi}}
\end{aligned}
$$

#### Application to convection equation

Let us consider the linear convection equation:

$$
    \pdv{u}{t} + a \pdv{u}{x} = 0 \implies \tilde{G} = e^{-\I a k \Delta t} \implies \veps_D = |G|, \veps_\phi = \frac{\Phi}{a k \Delta t} = \frac{a_\mrm{num}}{a}
$$



For the diffusion problem:
$$
    \pdv{u}{t} = \alpha \pmdv{u}{x}{2} \implies \tilde{G} = e^{-\alpha k^2}
$$