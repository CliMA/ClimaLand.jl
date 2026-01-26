# Multirate Runge Kutta

Given a problem with two components that operate at two rates:
```math
\frac{du}{dt} = f_F(u,t) + f_S(u,t)
```
where ``f_F`` is the _fast_ component, and ``f_S`` is the _slow_ component.

[SKAW2012](@cite) defines the following method.

Given an outer explicit Runge--Kutta scheme with tableau ``(a,b,c)``

We can define the stage values ``U^{(i)} = v_i(\tau_i)`` as the solution to the _inner_ ODE problem
```math
\frac{dv_i}{d\tau}
  = \sum_{j=1}^i  \frac{a_{ij} - a_{i-1,j}}{c_i - c_{i-1}}  f_S (U^{(j)}, \tau_j)
    + f_F(v_i, \tau),
\quad \tau \in [\tau_{i-1}, \tau]
```
where ``\tau_i = t + \Delta t c_i``, with initial condition ``v_i(\tau_{i-1}) = U^{(i-1)}``. If ``c_i == c_{i-1}``, we can treat it as a correction step:
```math
U^{(i)} = U^{(i-1)} + \Delta t \frac{\sum_{j=1}^i (a_{ij} - a_{i-1,j})}{c_i - c_{i-1}} f_S (U^{(j)}, \tau_i)
```
The final summation stage treating analogously to ``i=N+1``, with ``a_{N+1,j} = b_j`` and ``c_{N+1} = 1``.

## Low-storage

If using a low-storage Runge--Kutta method is used as an outer solver, then this reduces to
```math
\frac{dv_i}{d\tau}
  =  \frac{B_{i-1}}{c_i - c_{i-1}} dU_S^{(i-1)}
    + f_F(v_i, \tau),
\quad \tau \in [\tau_{i-1}, \tau]
```
where
```math
dU_S^{(i)} = f_S(U^{(i)}, \tau_i) + A_i dU_S^{(i-1)}
```

## Multirate Infinitesimal Step (MIS)

Multirate Infinitesimal Step (MIS) methods ([WKG2009](@cite), [KW2014](@cite))

```math
\begin{aligned}
v_i (0)
  &= u^n + \sum_{j=1}^{i-1} \alpha_{ij} (U^{(j)} - u^n)
\\
\frac{dv_i}{d\tau}
  &= \sum_{j=1}^{i-1} \frac{\gamma_{ij}}{d_i \Delta t} (U^{(j)} - u^n)
    + \sum_{j=1}^i \frac{\beta_{ij}}{d_i} f_S (U^{(j)}, t + \Delta t c_i)
    + f_F(v_i, t^n +  \Delta t \tilde c_i + \frac{c_i - \tilde c_i}{d_i} \tau),
\quad \tau \in [0, \Delta t d_i]
\\
U^{(i)} &= v_i(\Delta t d_i)
\end{aligned}
```
The method is defined in terms of the lower-triangular matrices ``\alpha``,
``\beta`` and ``\gamma``, with ``d_i = \sum_j \beta_{ij}``,
``c_i = (I - \alpha - \gamma)^{-1} d`` and ``\tilde c = \alpha c``.


## Wicker Skamarock

[WS1998](@cite) and [WS2002](@cite) define RK2 and RK3 multirate schemes:
```math
\begin{aligned}
v_i (t) &= u^n
\\
\frac{dv_i}{d\tau}
  &= f_S (U^{(i-1)}, t + \Delta t c_{i-1})
    + f_F(v_i, \tau),
\quad \tau \in [t, t+ \Delta t c_i ]
\\
U^{(i)} &= v_i(t + \Delta t c_i)
\end{aligned}
```
which corresponds to an MIS method with ``\alpha = \gamma = 0`` and ``\beta = \operatorname{diag}(c)``.