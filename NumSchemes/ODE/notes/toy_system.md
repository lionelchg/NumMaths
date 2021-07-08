# Toy differential system

$$
\frac{dY}{dt} = f(t, Y) = A Y \qquad \frac{dW}{dt} = D W \qquad Y = PW \ D = P^{-1}A P
$$

Let us take the following form for $D$ with $\alpha > 1$:

$$
D = \begin{pmatrix}
1 & \cdots & \cdots & \cdots \\
0 & \alpha & \cdots & \cdots \\
0 & 0 & \alpha^2 & \cdots \\
0 & 0 & 0 & - 1 \\
0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 \\
\end{pmatrix}
$$