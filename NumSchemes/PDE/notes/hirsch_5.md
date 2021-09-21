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

To be able to analyze a scheme associted to specific boundary conditions the matrix $S$ needs to be diagonalized. It is a complicated task, either tedious to do by hand or very costly numerically, so that periodic boundary conditions can be simplified so that the eigenvectors are the Fourier modes

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