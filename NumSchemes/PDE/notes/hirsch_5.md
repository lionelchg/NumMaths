# Numerical Computation of Internal and External Flows notes

## Chapter 9: Time Integration Methods for Space-discretized Equations

### Matrix representation of transport

In numerical simulations of flows we are interested in laws of the form

$$
\pdv{u}{t} + \nabla \cdot \vb{F} = 0
$$

Having a spatial discretization scheme expressed in matrix form with boundary conditions yields the following system of equations

$$
\dv{U}{t} = SU + Q
$$

The matrix $S$ which contains the discretization scheme is supposed to be diagonalisable with $\Omega_j, j \in \intint{1}{N}$ its eigenvalues so that

$$
\Omega = T^{-1} S T
$$

where $T$ is the transition matrix. The eigenvectors of the matrix, denoted $V^{(j)}(\vb{x})$ form a basis of the discretized vector space on which the values are defined. Hence we have the modal decomposition

$$
\bar{U}(t, \vb{x}) = \sum_{j=1}^N \bar{U}_j(t) V^{(j)}(\vb{x})
$$

After expanding the non-homogeneous term in the eigenvectors base $Q = \sum_j Q_j V^{(j)}$, each mode follows the following differential equation

$$
\dv{\bar{U}_j}{t} = \Omega_j \bar{U}_j(t) + Q_j
$$

The solution of this differential equation is simply

$$
\bar{U}_j(t) = U_j^0 e^{\Omega_j t} + \frac{Q_j}{\Omega_j} (e^{\Omega_j t} - 1)
$$

### Stability condition

The solution of the system is thus a sum over all the modal solutions. The system of ODEs is well-posed or stable if it is bounded as time advances. For this solution not to explode in time the real part of the eigenvalues need to be negative hence

$$
\Re(\Omega_j) \leq 0 \quad \forall j \in \intint{1}{N}
$$

and if an eigenvalue is zero, it has to be a simple eigenvalue.

### Matrix method and Fourier modes

To be able to analyze a scheme associted to specific boundary conditions the matrix $S$ needs to be diagonalized. It is a complicated task, either tedious to do by hand or very costly numerically, so that periodic boundary conditions can be assumed so that the eigenvectors are the Fourier modes

$$
V^{(j)}(x) = e^{\I k_j x}
$$

Inserting these eigenvectors back in the operator matrix we can find the eigenvalues associated to a specified scheme

$$
S e^{\I i \phi_j} = \Omega(\phi_j) e^{\I i \phi_j}
$$

### Amplification factor

By assuming $Q = 0$ the exact solution of the semi-discretized scheme is

$$
\bar{U}(n\Delta t) = \sum_{j=1}^N \bar{U}_{Tj}(n\Delta t) V^{(j)} = \sum_{j=1}^N U_j^0 e^{\Omega_j n \Delta t} V^{(j)}
$$

The amplitudes of each of the eigenvectors evole in time such that

$$
\bar{U}_{Tj}(n\Delta t) = e^{\Omega_j \Delta t} \bar{U}_{Tj}((n-1)\Delta t)
$$

So that the amplification factor of the semi-discretized scheme can be defined as

$$
G(\Omega) = e^{\Omega \Delta t}
$$

Hence a single mode can be isolated to study the behavior of the scheme in time. This leads to the canonical form of the modal equation

$$
\dv{w}{t} = \Omega w
$$

### Analysis of Time Integration schemes

We introduce the time advancement operator as

$$
E w^n = w^{n+1}
$$

The time integration scheme is defined generally as

$$
w^{n+1} = P(E, \Omega \Delta t) w^n
$$

where the operator $P$ is the numerical amplification factor of the time integration, depending only on $\Omega \Delta t$. Stability requires that this time integration amplification factor is bounded when $n \rightarrow \infty$ and $\Delta t \rightarrow 0$ for $n\Delta t$ fixed. For a selected norm and some specific $T$ we should have

$$
||P^n|| < K \quad \mrm{for} \: 0 < n \Delta t < T
$$

The norm of $P$ is difficult to analyze and a necessary but not always sufficient condition is rather applied from a local mode analysis. The eigenvalues of the operator $P$ are defined as

$$
z_P = P(z_P, \Omega \Delta t)
$$

and the condition is simply

$$
|z_P| \leq 1
$$

for all eigenvalues $\Omega$. $z_P$ is a function of $\Omega \Delta t$ so the condition of stability must be fulfilled as $\Omega \Delta t$ covers its spectrum. When the number of solutions is greater than one, there is one root, called the principal solution, which must tend to 1 when the time step tends to zero and which represents an approximation to the physical behavior. The other spurious solutions, which represent non-physical behavior, can affect the stability of the solution although being non-physical.

### Error analysis of time and space discretizations

Replacing $P$ by $z_P$ local modes evolve in time according to

$$
w^n = z_P^n w^0
$$

This amplification factor of the time integration scheme $z_P$, associated to the eigenvalue $\Omega$, must be an approximation of the exact amplification factor of the semi-discretized scheme $\bar{G} = \exp(\Omega \Delta t)$ so that

$$
z_P(\Omega \Delta t) \simeq  \exp(\Omega \Delta t) = 1 + \Omega \Delta t + \frac{\Omega \Delta t}{2!} + \frac{\Omega \Delta t}{3!} + \ldots
$$

Hence when writing the Taylor expansion of $z_P$, the first time that is different from the Taylor expansion of the exponential defines the order of the time integration scheme.

The time integration error is defined as

$$
\veps^T = \frac{z_P}{\bar{G}} = \frac{|z_P|}{\exp(\mrm{Re}(\Omega \Delta t))} \frac{\exp(\I \Phi_p)}{\exp(\mrm{Re}(\Omega \Delta t))}
$$

Likewise the space integration error is

$$
\veps^S = \frac{\bar{G}}{\tilde{G}}
$$

where $\tilde{G}$ is the exact amplification factor of the initial differential equation.

### Runge-Kutta schemes

The low-storage Runge-Kutta (RK) schemes are defined as

$$
\begin{aligned}
&U^{(1)} = U^n \\
&\cdots \\
&U^{(j)} = U^n + \alpha_j H^{(j-1)}\\
&\cdots \\
&U^{(K)} = U^n + \alpha_K H^{(K-1)}\\
&U^{n+1} = U^n + \Delta t \sum_{k=1}^K \beta_k H^{(k)}
\end{aligned}
$$

where $\sum_k \beta_k = 1$ for consistency and each $U^{(j)}$ is called a RK stage and $K$ the number of stages. In AVIP the choice of $\beta_K = 1$, $\beta_j = 0 \quad j=1 \ldots k-1$ is made. The RK1 to RK4 methods are defined as

| Time integration method | $\alpha_2$ | $\alpha_3$ | $\alpha_4$ | $\beta_K$ |
| ----------------------- | ---------- | ---------- | ---------- | --------- |
| RK1                     | 0          | 0          | 0          | 1         |
| RK2                     | 1/2        | 0          | 0          | 1         |
| RK3                     | 1/3        | 1/2        | 0          | 1         |
| RK4                     | 1/4        | 1/3        | 1/2        | 1         |