# Low-storage Runge--Kutta methods

LSRK methods are self-starting, with ``U^{(1)} = u^n``, and then using stage updates of the form
```math
\begin{aligned}
dU^{(i)} &= f(U^{(i)}, t + c_i \Delta t) + A_i dU^{(i-1)}\\
U^{(i+1)} &= U^{(i)} + \Delta t B_i dU^{(i)}
\end{aligned}
```
where ``A_1 = c_1 = 0`` (implying ``dU^{(1)} = f(u^n, t)``), with the value at the next step being the ``N+1``th stage value ``u^{n+1} = U^{(N+1)})``.

This allows the updates to be performed with only two copies of the state vector (so long as `f` can be evaluated in incrementing form).

It can be written as an RK scheme with Butcher tableau coefficients defined by the recurrences
```math
\begin{aligned}
a_{i,i} &= 0 \\
a_{i,j} &= B_j + A_{j+1} a_{i,j+1}\\
b_N &= B_N \\
b_i &= B_i + A_{i+1} b_{i+1}
\end{aligned}
```
or equivalently
```math
\begin{aligned}
a_{j,j} &= 0 \\
a_{i,j} &= a_{i-1,j} + B_{i-1} \prod_{k=j+1}^{i-1} A_k
\end{aligned}
```
with ``b_j`` treated analogously as ``a_{N+1,j}``.
