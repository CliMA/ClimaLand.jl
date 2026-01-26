# ODE Solvers

## Standard IMEX ARK

An *ordinary differential equation* (ODE) is an equation of the form

```math
\frac{d}{dt}u(t) = T(u(t), t),
```

where ``u(t)`` is called the *state* at time ``t``, and ``T(u(t), t)`` is called the *tendency* of the state at time ``t``.

The simplest method for numerically solving this equation with a finite *timestep* ``\Delta t`` is the *forward Euler* method, in which the equation is approximated as

```math
\frac{u(t + \Delta t) - u(t)}{\Delta t} \approx T(u(t), t).
```

Given the value of ``u_0 = u(t_0)`` at some time ``t_0``, the approximation implies that ``u(t_0 + \Delta t) \approx \hat{u}``, where

```math
\hat{u} = u_0 + \Delta t T(u_0, t_0).
```

An alternative approximation is given by the *backward Euler* method:

```math
\frac{u(t + \Delta t) - u(t)}{\Delta t} \approx T(u(t + \Delta t), t + \Delta t).
```

With this approximation, ``u(t_0 + \Delta t) \approx \hat{u}``, where ``\hat{u}`` is now the solution to the equation

```math
u_0 + \Delta t T(\hat{u}, t_0 + \Delta t) - \hat{u} = 0.
```

Unlike the forward Euler method, in which ``\hat{u}`` is directly computed based on the known state ``u_0``, the backward Euler method involves solving a root equation in order to obtain the value of ``\hat{u}``. So, the forward Euler method is called an *explicit* method, whereas the backward Euler method is called an *implicit* method.

In general, ``T`` can be a complicated nonlinear function of ``u(t)`` and ``t``, and it is usually not possible to solve the backward Euler method's implicit equation for ``\hat{u}`` analytically. Instead, it is often necessary to use an iterative root-finding algorithm like Newton's method to solve for ``\hat{u}``. Although this is more computationally expensive than using the forward Euler method to directly obtain ``\hat{u}`` from ``u_0``, it is often necessary to deal with the problem of *stiffness*.

Roughly speaking, a tendency ``T`` is stiff when the forward Euler method requires a relatively small timestep ``\Delta t`` to obtain a reasonably accurate solution, where "relatively small" means that it is smaller than one would expect based on the rate at which ``u(t)`` changes over time. When ``T`` is not stiff, a small timestep is only required when ``u(t)`` changes quickly with respect to ``t``, and a larger timestep can be used when ``u(t)`` changes slowly. The backward Euler method is more *stable* than the forward Euler method, which means that it is usually able to take larger timesteps, especially when ``T`` is stiff.

As a compromise between the simplicity of the forward Euler method and the stability of the backward Euler method, it is common to use the *implicit-explicit* (IMEX) Euler method. This involves splitting ``T`` into an explicit tendency ``T_{\text{exp}}`` and an implicit tendency ``T_{\text{imp}}``, so that ``T(u(t), t) = T_{\text{exp}}(u(t), t) + T_{\text{imp}}(u(t), t)``, and making the approximation

```math
\frac{u(t + \Delta t) - u(t)}{\Delta t} \approx T_{\text{exp}}(u(t), t) + T_{\text{imp}}(u(t + \Delta t), t + \Delta t).
```

With this approximation, the implicit equation for ``\hat{u}`` becomes

```math
u_0 + \Delta t T_{\text{exp}}(u_0, t_0) + \Delta t T_{\text{imp}}(\hat{u}, t_0 + \Delta t) - \hat{u} = 0.
```

Note that ``T_{\text{exp}}`` is evaluated explicitly at the known state ``u_0``, while ``T_{\text{imp}}`` is evaluated implicitly at the unknown state ``\hat{u}``.

Both the forward and backward Euler methods (and, by extension, the IMEX Euler method) are *first-order* methods, which means that the errors of their approximations are proportional to ``\Delta t`` when ``\Delta t`` is sufficiently close to 0.[^1] In order to achieve a reasonable accuracy with fewer timesteps, it is common to use *higher-order* methods, where a method of order ``p`` will have an error that is proportional to ``(\Delta t)^p`` for small values of ``\Delta t``. The simplest higher-order generalization of the forward and backward Euler methods is a *Runge-Kutta* method, in which there are ``s`` *stages* ``U_1, U_2, \ldots, U_s`` that satisfy

[^1]: \
    More precisely, the *local truncation error* of the forward and backward Euler methods after a single timestep is ``O\left((\Delta t)^2\right)``, which means that ``|u(t_0 + \Delta t) - \hat{u}| < C (\Delta t)^2`` for all ``\Delta t < D``, where ``C`` and ``D`` are some constants. On the other hand, the *global truncation error* after taking enough timesteps to go from ``t_0`` to some ``t_1 > t_0`` is ``O(\Delta t)``. This is because there are ``(t_1 - t_0) / \Delta t`` timesteps between ``t_0`` and ``t_1``, and, if each timestep has a local truncation error of ``O\left((\Delta t)^2\right)``, then the error after ``O\left((\Delta t)^{-1}\right)`` timesteps must be ``O(\Delta t)``. In general, for a Runge-Kutta method (or ARK method) of order ``p``, the local truncation error is ``O\left((\Delta t)^{p + 1}\right)``, and the global truncation error is ``O\left((\Delta t)^{p\vphantom{1}}\right)``.

```math
U_i = u_0 + \Delta t \sum_{j = 1}^s a_{i,j} T(U_j, t_0 + \Delta t c_j),
```

and ``u(t_0 + \Delta t)`` is approximated as

```math
\hat{u} = u_0 + \Delta t \sum_{i = 1}^s b_i T(U_i, t_0 + \Delta t c_i).
```

The coefficients ``a_{i,j}``, ``b_i``, and ``c_i`` of a Runge-Kutta method can be summarized in a *Butcher tableau*:

```math
\begin{array}{c|c c c c c} c_1 & a_{1,1} & a_{1,2} & \cdots & a_{1,s - 1} & a_{1,s} \\ c_2 & a_{2,1} & a_{2,2} & \cdots & a_{2,s - 1} & a_{2,s} \\ c_3 & a_{3,1} & a_{3,2} & \cdots & a_{3,s - 1} & a_{3,s} \\ \vdots & \vdots & \vdots & \ddots & \vdots & \vdots \\ c_s & a_{s,1} & a_{s,2} & \cdots & a_{s,s - 1} & a_{s,s} \\ \hline & b_1 & b_2 & \cdots & b_{s - 1} & b_s \end{array}
```

Since solving a system of ``s`` coupled equations for the stages ``U_i`` is usually impractical, many of the ``a_{i,j}`` coefficients are set to 0 in commonly used Runge-Kutta methods. When ``a_{i,j} = 0`` for all ``j \geq i``, the equation for ``U_i`` simplifies to the explicit formula

```math
U_i = u_0 + \Delta t \sum_{j = 1}^{i - 1} a_{i,j} T(U_j, t_0 + \Delta t c_j).
```

This is called an explicit Runge-Kutta (ERK) method. When ``a_{i,j} = 0`` for all ``j > i``, the implicit equation for ``U_i`` becomes uncoupled from the equations for the other stages, and it can be rewritten as

```math
u_0 + \Delta t \sum_{j = 1}^{i - 1} a_{i,j} T(U_j, t_0 + \Delta t c_j) + \Delta t a_{i,i} T(U_i, t_0 + \Delta t c_i) - U_i = 0.
```

Since the only unknown tendency in this equation comes from the diagonal coefficient ``a_{i,i}``, this is called a *diagonally implicit* Runge-Kutta (DIRK) method. There are many different categories of DIRK methods, such as *explicit first stage* DIRK (EDIRK), where ``a_{1,1} = 0``, *singly* DIRK (SDIRK), where ``a_{1,1} = a_{2,2} = \ldots = a_{s,s} = \gamma`` for some constant ``\gamma``, and ESDIRK, where ``a_{1,1} = 0`` and ``a_{2,2} = \ldots = a_{s,s} = \gamma``.

Just as the forward and backward Euler methods can be combined into the IMEX Euler method, two Runge-Kutta methods can be combined into an *additive* Runge-Kutta (ARK) method. If the first method is explicit and the second is implicit, the result is an IMEX ARK method. We will only be considering IMEX ARK methods where the implicit part is DIRK. If the DIRK method has coefficients ``a_{i,j}``, ``b_i``, and ``c_i``, and the ERK method has coefficients ``\tilde{a}_{i,j}``, ``\tilde{b}_i``, and ``\tilde{c}_i``, the IMEX ARK method approximates ``u(t_0 + \Delta t)`` as

```math
\hat{u} = u_0 + \Delta t \sum_{i = 1}^s \left(\tilde{b}_i T_{\text{exp}}(U_i, t_0 + \Delta t \tilde{c}_i) + b_i T_{\text{imp}}(U_i, t_0 + \Delta t c_i)\right),
```

where the implicit equation for ``U_i`` is now

```math
u_0 + \Delta t \sum_{j = 1}^{i - 1} \left(\tilde{a}_{i,j} T_{\text{exp}}(U_j, t_0 + \Delta t \tilde{c}_j) + a_{i,j} T_{\text{imp}}(U_j, t_0 + \Delta t c_j)\right) + \Delta t a_{i,i} T_{\text{imp}}(U_i, t_0 + \Delta t c_i) - U_i = 0.
```

## Adding DSS

It is often necessary to filter the state so that it satisfies some particular constraint before it is used to evaluate a tendency. In our case, the state is a collection of values defined across a spatially discretized domain. When we use a *continuous Galerkin* (CG) spectral element discretization, we must ensure that the state is continuous across element boundaries before we use it to compute any tendency. We can turn any state that is discontinuous across element boundaries into a continuous one by applying a *direct stiffness summation* (DSS) to it. Applying DSS to ``\hat{u}`` in the IMEX ARK method is straightforward:

```math
\hat{u} = \textrm{DSS}\left(u_0 + \Delta t \sum_{i = 1}^s \left(\tilde{b}_i T_{\text{exp}}(U_i, t_0 + \Delta t \tilde{c}_i) + b_i T_{\text{imp}}(U_i, t_0 + \Delta t c_i)\right)\right).
```

Applying DSS to each stage ``U_i`` is a bit trickier. Ideally, we would use the equation

```math
\textrm{DSS}\left(\begin{aligned} & u_0 + \Delta t \sum_{j = 1}^{i - 1} \left(\tilde{a}_{i,j} T_{\text{exp}}(U_j, t_0 + \Delta t \tilde{c}_j) + a_{i,j} T_{\text{imp}}(U_j, t_0 + \Delta t c_j)\right) + \\ & \Delta t a_{i,i} T_{\text{imp}}(U_i, t_0 + \Delta t c_i) \end{aligned}\right) - U_i = 0,
```

since this would ensure that the implicit tendency ``T_{\text{imp}}`` gets evaluated at a continuous stage ``U_i``. However, this equation is more challenging to solve than the original one, since it involves applying DSS, which is usually a more complicated function than ``T_{\text{imp}}``, to an unknown quantity. So, we instead use the equation

```math
\begin{aligned} \textrm{DSS}\left(u_0 + \Delta t \sum_{j = 1}^{i - 1} \left(\tilde{a}_{i,j} T_{\text{exp}}(U_j, t_0 + \Delta t \tilde{c}_j) + a_{i,j} T_{\text{imp}}(U_j, t_0 + \Delta t c_j)\right)\right) + & \\ \Delta t a_{i,i} T_{\text{imp}}(U_i, t_0 + \Delta t c_i) - U_i &= 0. \end{aligned}
```

This equation is identical to the previous one only when ``T_{\text{imp}}`` preserves the effects of DSS. That is to say, if evaluating ``T_{\text{imp}}`` at a continuous stage ``U_i`` produces a continuous tendency, then it is not necessary to apply DSS after adding this tendency to a continuous quantity. This also assumes that the root-finding algorithm used to compute ``U_i`` preserves the effects of DSS, which will usually be the case as long as ``T_{\text{imp}}`` does so. In general, one could replace DSS with a filter that enforces some other constraint, as long as that constraint is preserved by ``T_{\text{imp}}``. The constraint must also be preserved by addition and by multiplication with a constant, which are the only other operations used to compute ``U_i`` after applying the filter.

Note that this method of enforcing a constraint is not mathematically rigorous, as it does not necessarily maintain the convergence properties of the IMEX ARK method. This is because DSS effectively acts like a tendency on each stage ``U_i``, but, unlike the actual tendencies ``T_{\text{exp}}`` and ``T_{\text{imp}}``, the effects of DSS on the previous stages ``U_1, U_2, \ldots, U_{i - 1}`` are not accounted for when computing ``U_i``. The proper way to enforce a constraint would be to extend the ODE to a *differential-algebraic equation* (DAE) by adding an algebraic equation ``F(u(t), t) = 0`` that can only be satisfied when ``u(t)`` obeys the constraint, and then solving the DAE using an appropriate numerical method. However, this route would be significantly more computationally expensive than our method. Moreover, we have observed that the effects of DSS are sufficiently small to not noticeably disrupt the convergence of IMEX ARK methods in our test cases, though this will not necessarily be the case if DSS is replaced with another filter.

## Adding a Limiter

In addition to filtering the state, it is often necessary to filter part of the tendency, but in a way that depends on the state to which the tendency is being added. In our case, because we use a spectral element discretization, part of the explicit tendency can add spurious oscillations to each stage. We can limit these oscillations by using a *monotonicity-preserving limiter* (see [GTS2014](@cite)), which is a function ``\textrm{lim}_{u(t)}`` that, when used with an appropriate tendency ``T_{\text{lim}}(u(t), t)`` and a sufficiently small constant ``C``, satisfies the inequalities

```math
\begin{aligned} \textrm{min}\biggl(\textrm{lim}_{u(t)}\bigl(u(t) + C T_{\text{lim}}(u(t), t)\bigr)\biggr) &\geq \textrm{min}\bigl(u(t)\bigr) \text{ and} \\ \textrm{max}\biggl(\textrm{lim}_{u(t)}\bigl(u(t) + C T_{\text{lim}}(u(t), t)\bigr)\biggr) &\leq \textrm{max}\bigl(u(t)\bigr). \end{aligned}
```

In other words, applying the limiter to a state incremented by the limited tendency, ``u(t) + C T_{\text{lim}}(u(t), t)``, ensures that the extrema of the incremented state do not exceed the extrema of the unincremented state, ``u(t)``. Note that the process of incrementing a state in this way is used by the forward Euler method with ``\Delta t = C`` and ``T = T_{\text{lim}}``, which is why it is called an *Euler step*. Since spurious oscillations usually cause the extrema of the state to grow, this is an effective mechanism for eliminating such oscillations. Note that the limiter can only preserve monotonicity by comparing the incremented state to the unincremented state, which means that the effects of the limiter are a function of the unincremented state (hence the subscript in ``\textrm{lim}_{u(t)}``).

Unfortunately, there is no mathematically correct way to incorporate the use of a limiter into a general IMEX ARK method. The most straightforward approach is to split ``T_{\text{exp}}`` into two tendencies ``T_{\text{exp}}`` and ``T_{\text{lim}}``, to modify the equation for ``\hat{u}`` to

```math
\begin{aligned} \hat{u} ={} & \textrm{lim}_{u_0}\left(u_0 + \Delta t \sum_{i = 1}^s \tilde{b}_i T_{\text{lim}}(U_i, t_0 + \Delta t \tilde{c}_i)\right) + \\ & \Delta t \sum_{i = 1}^s \left(\tilde{b}_i T_{\text{exp}}(U_i, t_0 + \Delta t \tilde{c}_i) + b_i T_{\text{imp}}(U_i, t_0 + \Delta t c_i)\right), \end{aligned}
```

and to modify the equation for ``U_i`` to

```math
\begin{aligned} \textrm{lim}_{u_0}\left(u_0 + \Delta t \sum_{j = 1}^{i - 1} \tilde{a}_{i,j} T_{\text{lim}}(U_j, t_0 + \Delta t \tilde{c}_j)\right) + & \\ \Delta t \sum_{j = 1}^{i - 1} \left(\tilde{a}_{i,j} T_{\text{exp}}(U_j, t_0 + \Delta t \tilde{c}_j) + a_{i,j} T_{\text{imp}}(U_j, t_0 + \Delta t c_j)\right) + & \\ \Delta t a_{i,i} T_{\text{imp}}(U_i, t_0 + \Delta t c_i) - U_i &= 0. \end{aligned}
```

Not only does this approach not maintain the convergence properties of the IMEX ARK method (for the same reason as DSS above),[^2] but it also does not even use the limiter in a way that allows it to preserve monotonicity. That is, the argument of ``\textrm{lim}_{u_0}`` does not have the form needed to satisfy the min and max constraints given above, since it involves evaluating the limited tendency at states that are not the unincremented state. As luck would have it, we have observed that the limiter, even when used in such an incorrect manner, is still able to substantially reduce spurious oscillations in our test cases, though it is not able to perform as well as it would if it were used correctly.

[^2]: \
    By [Godunov's theorem](https://en.wikipedia.org/wiki/Godunov%27s_theorem), no monotonicity-preserving linear numerical method can have an order greater than 1. Since ``\textrm{lim}_{u_0}`` will usually be a nonlinear function, this is not a linear numerical method. However, it is a rough approximation of the unmodified ARK method, so it is likely that Godunov's theorem will still apply; i.e., we do not expect to observe an order greater than 1 when using a limiter.

In order to use the limiter "more correctly", we constrain the ERK Butcher tableau coefficients to have the following form:

```math
\begin{array}{c|c c c c c} \tilde{c}_1 & 0 & 0 & \cdots & 0 & 0 \\ \tilde{c}_2 & \beta_1 & 0 & \cdots & 0 & 0 \\ \tilde{c}_3 & \beta_1 \beta_2 & \beta_2 & \cdots & 0 & 0 \\ \vdots & \vdots & \vdots & \ddots & \vdots & \vdots \\ \tilde{c}_s & \beta_1 \beta_2 \ldots \beta_{s - 1} & \beta_2 \ldots \beta_{s - 1} & \cdots & \beta_{s - 1} & 0 \\ \hline & \beta_1 \beta_2 \ldots \beta_{s - 1} \beta_s & \beta_2 \ldots \beta_{s - 1} \beta_s & \cdots & \beta_{s - 1} \beta_s & \beta_s \end{array}
```

If all of the coefficients ``\beta_1, \beta_2, \ldots, \beta_s`` are nonnegative, and if the resulting ERK method converges, it is called a *strong stability preserving* Runge-Kutta (SSPRK) method.[^3] Moreover, if the ERK method is SSPRK, then the overall IMEX ARK method is called an IMEX SSPRK method.

[^3]: \
    This is really a *low-storage* SSPRK method that only requires two registers for storing states, using one of the registers to store ``u_0``; this is an example of a ``2N^*`` method (see [HR2018arxiv](@cite)). A general SSPRK method has just as many independent coefficients as a general ERK method, though its coefficients are used somewhat differently, in a way that makes them more amenable to limiters. We restrict ourselves to ``s`` independent coefficients ``\beta_i`` merely for the sake of simplicity. In the future, we might want to generalize to arbitrary SSPRK methods.

Now, in order to simplify our notation, we will define ``s + 1`` values ``\tilde{U}_1, \tilde{U}_2, \ldots, \tilde{U}_{s + 1}``, where

```math
\tilde{U}_i = \begin{cases} \displaystyle u_0 + \Delta t \sum_{j = 1}^{i - 1} \tilde{a}_{i,j} T_{\text{exp}}(U_j, t_0 + \Delta t \tilde{c}_j) & i < s + 1 \\ \displaystyle u_0 + \Delta t \sum_{i = 1}^s \tilde{b}_i T_{\text{exp}}(U_i, t_0 + \Delta t \tilde{c}_i) & i = s + 1 \end{cases}.
```

This allows us to rewrite the IMEX ARK equation for ``\hat{u}`` as

```math
\hat{u} = \tilde{U}_{s + 1} + \Delta t \sum_{i = 1}^s b_i T_{\text{imp}}(U_i, t_0 + \Delta t c_i),
```

and to rewrite the equation for ``U_i`` as

```math
\tilde{U}_i + \Delta t \sum_{j = 1}^{i - 1} a_{i,j} T_{\text{imp}}(U_j, t_0 + \Delta t c_j) + \Delta t a_{i,i} T_{\text{imp}}(U_i, t_0 + \Delta t c_i) - U_i = 0.
```

If we constrain the IMEX ARK method to an IMEX SSPRK method, we can express ``\tilde{b}_i`` as

```math
\tilde{b}_i = \begin{cases} \beta_s \tilde{a}_{s,i} & i < s \\ \beta_s & i = s \end{cases},
```

which means that

```math
\begin{aligned} \tilde{U}_{s + 1} ={} & u_0 + \Delta t \sum_{i = 1}^{s - 1} \tilde{b}_i T_{\text{exp}}(U_i, t_0 + \Delta t \tilde{c}_i) + \Delta t \tilde{b}_s T_{\text{exp}}(U_s, t_0 + \Delta t \tilde{c}_s) = \\ & u_0 + \Delta t \beta_s \sum_{i = 1}^{s - 1} \tilde{a}_{s,i} T_{\text{exp}}(U_i, t_0 + \Delta t \tilde{c}_i) + \Delta t \beta_s T_{\text{exp}}(U_s, t_0 + \Delta t \tilde{c}_s) = \\ & (1 - \beta_s) u_0 + \beta_s \left(\tilde{U}_s + \Delta t T_{\text{exp}}(U_s, t_0 + \Delta t \tilde{c}_s)\right). \end{aligned}
```

In addition, for all ``i > 1``,

```math
\tilde{a}_{i,j} = \begin{cases} \beta_{i - 1} \tilde{a}_{i - 1,j} & j < i - 1 \\ \beta_{i - 1} & j = i - 1 \end{cases},
```

which means that

```math
\begin{aligned} \tilde{U}_i ={} & u_0 + \Delta t \sum_{j = 1}^{i - 2} \tilde{a}_{i,j} T_{\text{exp}}(U_j, t_0 + \Delta t \tilde{c}_j) + \Delta t \tilde{a}_{i,i - 1} T_{\text{exp}}(U_{i - 1}, t_0 + \Delta t \tilde{c}_{i - 1}) = \\ & u_0 + \Delta t \beta_{i - 1} \sum_{j = 1}^{i - 2} \tilde{a}_{i - 1,j} T_{\text{exp}}(U_j, t_0 + \Delta t \tilde{c}_j) + \Delta t \beta_{i - 1} T_{\text{exp}}(U_{i - 1}, t_0 + \Delta t \tilde{c}_{i - 1}) = \\ & (1 - \beta_{i - 1}) u_0 + \beta_{i - 1} \left(\tilde{U}_{i - 1} + \Delta t T_{\text{exp}}(U_{i - 1}, t_0 + \Delta t \tilde{c}_{i - 1})\right). \end{aligned}
```

Since ``\tilde{U}_1 = u_0``, constraining the IMEX ARK method to an IMEX SSPRK method allows us to express ``\tilde{U}_i`` as

```math
\tilde{U}_i = \begin{cases} u_0 & i = 1 \\ (1 - \beta_{i - 1}) u_0 + \beta_{i - 1} \left(\tilde{U}_{i - 1} + \Delta t T_{\text{exp}}(U_{i - 1}, t_0 + \Delta t \tilde{c}_{i - 1})\right) & i > 1 \end{cases}.
```

To incorporate the use of a limiter into the IMEX SSPRK method, we split ``T_{\text{exp}}`` into ``T_{\text{exp}}`` and ``T_{\text{lim}}``, and we modify the equation for ``\tilde{U}_i`` to

```math
\tilde{U}_i = \begin{cases} u_0 & i = 1 \\ \begin{aligned} & (1 - \beta_{i - 1}) u_0 +{} \\ & \quad\beta_{i - 1} \left(\textrm{lim}_{\tilde{U}_{i - 1}}\left(\tilde{U}_{i - 1} + \Delta t T_{\text{lim}}(U_{i - 1}, t_0 + \Delta t \tilde{c}_{i - 1})\right) + \Delta t T_{\text{exp}}(U_{i - 1}, t_0 + \Delta t \tilde{c}_{i - 1})\right) \end{aligned} & i > 1 \end{cases}.
```

In this equation, the limiter is being applied to a limited tendency evaluated at ``U_{i - 1}``, but with ``\tilde{U}_{i - 1}`` as the unincremented state. If there is no implicit tendency, so that ``T_{\text{imp}}`` is always 0, then ``U_{i - 1} = \tilde{U}_{i - 1}``, and the limiter is able to properly preserve monotonicity. On the other hand, if there is an implicit tendency, then the limiter will not necessarily preserve monotonicity. That is, the limiter is guaranteed to function properly when the limited tendency is used in a sequence of Euler steps.

## Summary

We will now summarize our IMEX methods when using both DSS and a limiter.

#### IMEX ARK

Our general IMEX ARK method is defined by two Butcher tableaus:

```math
\begin{array}{c|c c c c c} \tilde{c}_1 & 0 & 0 & \cdots & 0 & 0 \\ \tilde{c}_2 & \tilde{a}_{2,1} & 0 & \cdots & 0 & 0 \\ \tilde{c}_3 & \tilde{a}_{3,1} & \tilde{a}_{3,2} & \cdots & 0 & 0 \\ \vdots & \vdots & \vdots & \ddots & \vdots & \vdots \\ \tilde{c}_s & \tilde{a}_{s,1} & \tilde{a}_{s,2} & \cdots & \tilde{a}_{s,s - 1} & 0 \\ \hline & \tilde{b}_1 & \tilde{b}_2 & \cdots & \tilde{b}_{s - 1} & \tilde{b}_s \end{array} \text{ and } \begin{array}{c|c c c c c} c_1 & a_{1,1} & 0 & \cdots & 0 & 0 \\ c_2 & a_{2,1} & a_{2,2} & \cdots & 0 & 0 \\ c_3 & a_{3,1} & a_{3,2} & \cdots & 0 & 0 \\ \vdots & \vdots & \vdots & \ddots & \vdots & \vdots \\ c_s & a_{s,1} & a_{s,2} & \cdots & a_{s,s - 1} & a_{s,s} \\ \hline & b_1 & b_2 & \cdots & b_{s - 1} & b_s \end{array}
```

Given ``u_0 = u(t_0)``, it approximates the value of ``u(t_0 + \Delta t)`` as

```math
\hat{u} = \textrm{DSS}\left(\begin{aligned} & \textrm{lim}_{u_0}\left(u_0 + \Delta t \sum_{i = 1}^s \tilde{b}_i T_{\text{lim}}(U_i, t_0 + \Delta t \tilde{c}_i)\right) + \\ & \Delta t \sum_{i = 1}^s \left(\tilde{b}_i T_{\text{exp}}(U_i, t_0 + \Delta t \tilde{c}_i) + b_i T_{\text{imp}}(U_i, t_0 + \Delta t c_i)\right) \end{aligned}\right),
```

where ``U_i`` is the solution to the equation

```math
\textrm{DSS}\left(\begin{aligned} & \textrm{lim}_{u_0}\left(u_0 + \Delta t \sum_{j = 1}^{i - 1} \tilde{a}_{i,j} T_{\text{lim}}(U_j, t_0 + \Delta t \tilde{c}_j)\right) + \\ & \Delta t \sum_{j = 1}^{i - 1} \left(\tilde{a}_{i,j} T_{\text{exp}}(U_j, t_0 + \Delta t \tilde{c}_j) + a_{i,j} T_{\text{imp}}(U_j, t_0 + \Delta t c_j)\right) \end{aligned}\right) + \Delta t a_{i,i} T_{\text{imp}}(U_i, t_0 + \Delta t c_i) - U_i = 0.
```

#### IMEX SSPRK

Our IMEX SSPRK method is defined by two Butcher tableaus:

```math
\begin{array}{c|c c c c c} \tilde{c}_1 & 0 & 0 & \cdots & 0 & 0 \\ \tilde{c}_2 & \beta_1 & 0 & \cdots & 0 & 0 \\ \tilde{c}_3 & \beta_1 \beta_2 & \beta_2 & \cdots & 0 & 0 \\ \vdots & \vdots & \vdots & \ddots & \vdots & \vdots \\ \tilde{c}_s & \beta_1 \beta_2 \ldots \beta_{s - 1} & \beta_2 \ldots \beta_{s - 1} & \cdots & \beta_{s - 1} & 0 \\ \hline & \beta_1 \beta_2 \ldots \beta_{s - 1} \beta_s & \beta_2 \ldots \beta_{s - 1} \beta_s & \cdots & \beta_{s - 1} \beta_s & \beta_s \end{array} \text{ and } \begin{array}{c|c c c c c} c_1 & a_{1,1} & 0 & \cdots & 0 & 0 \\ c_2 & a_{2,1} & a_{2,2} & \cdots & 0 & 0 \\ c_3 & a_{3,1} & a_{3,2} & \cdots & 0 & 0 \\ \vdots & \vdots & \vdots & \ddots & \vdots & \vdots \\ c_s & a_{s,1} & a_{s,2} & \cdots & a_{s,s - 1} & a_{s,s} \\ \hline & b_1 & b_2 & \cdots & b_{s - 1} & b_s \end{array}
```

Given ``u_0 = u(t_0)``, it approximates the value of ``u(t_0 + \Delta t)`` as

```math
\hat{u} = \textrm{DSS}\left(\tilde{U}_{s + 1} + \Delta t \sum_{i = 1}^s b_i T_{\text{imp}}(U_i, t_0 + \Delta t c_i)\right),
```

where ``U_i`` is the solution to the equation

```math
\textrm{DSS}\left(\tilde{U}_i + \Delta t \sum_{j = 1}^{i - 1} a_{i,j} T_{\text{imp}}(U_j, t_0 + \Delta t c_j)\right) + \Delta t a_{i,i} T_{\text{imp}}(U_i, t_0 + \Delta t c_i) - U_i = 0,
```

and where

```math
\tilde{U}_i = \begin{cases} u_0 & i = 1 \\ \begin{aligned} & (1 - \beta_{i - 1}) u_0 +{} \\ & \quad\beta_{i - 1} \left(\textrm{lim}_{\tilde{U}_{i - 1}}\left(\tilde{U}_{i - 1} + \Delta t T_{\text{lim}}(U_{i - 1}, t_0 + \Delta t \tilde{c}_{i - 1})\right) + \Delta t T_{\text{exp}}(U_{i - 1}, t_0 + \Delta t \tilde{c}_{i - 1})\right) \end{aligned} & i > 1 \end{cases}.
```

## Running Newton's Method

Every implicit equation for ``U_i`` has the form ``f_i(U_i) = 0``, where

```math
f_i(x) = R_i + \Delta t a_{i,i} T_{\text{imp}}(x, t_0 + \Delta t c_i) - x.
```

In this function, ``R_i``, ``\Delta t a_{i,i}``, and ``t_0 + \Delta t c_i`` are all quantities that do not depend on ``x``. The Jacobian of this function is

```math
W_i(x) = \frac{d}{dx}f_i(x) = \Delta t a_{i,i} J_{\text{imp}}(x, t_0 + \Delta t c_i) - 1,
```

where ``J_{\text{imp}}`` is the Jacobian of the implicit tendency,

```math
J_{\text{imp}}(x, t) = \frac{\partial}{\partial x}T_{\text{imp}}(x, t).
```

The value of ``U_i`` can be computed by running Newton's method with ``f = f_i`` and ``j = W_i``.

Note that "``W``" is used to denote the same quantity in [OrdinaryDiffEq.jl](https://github.com/SciML/OrdinaryDiffEq.jl/blob/v6.0.0/src/derivative_utils.jl).
