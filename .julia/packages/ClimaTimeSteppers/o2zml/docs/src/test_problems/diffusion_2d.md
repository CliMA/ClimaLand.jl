# 2D diffusion problem

TODO: cleanup

```
Due to conservation of (heat) energy, we know that
    c * ρ * ∂u/∂t = k * Δu + q,
where u is the temperature (K), c is the specific heat capacity (J/kg/K), ρ is
the density (kg/m^3), k is the thermal conductivity (W/m/K), and q is the rate
at which heat energy is added/removed (W/m^3).
We can simplify this PDE to
    ∂u/∂t = α * Δu + f,
where α = k/c/ρ is the thermal diffusivity (m^2/s) and f = q/c/ρ is the rate at
which heat energy is added/removed in units of temperature (K/s).
We will solve this PDE for u(x, y, t) over the domain (x, y) ∈ [0, l] × [0, l]
and t ≥ 0. For simplicity, we will use periodic boundary conditions:
    u(0, y, t) = u(l, y, t),
    u(x, 0, t) = u(x, l, t),
    ∇u(0, y, t) = ∇u(l, y, t), and
    ∇u(x, 0, t) = ∇u(x, l, t).
Also, for simplicity, we will assume that α is a constant.
Suppose that
    f = 0 and
    u(x, y, 0) = u₀(x, y).
The general solution to the PDE (obtained with separation of variables) is then
    u(x, y, t) =
        ∑_{n = 0}^∞ ∑_{m = 0}^∞ exp(-λₙₘ * t) * (
            φᶜᶜₙₘ(x, y) * ⟨φᶜᶜₙₘ(x, y), u₀(x, y)⟩ / ⟨φᶜᶜₙₘ(x, y), φᶜᶜₙₘ(x, y)⟩ +
            φᶜˢₙₘ(x, y) * ⟨φᶜˢₙₘ(x, y), u₀(x, y)⟩ / ⟨φᶜˢₙₘ(x, y), φᶜˢₙₘ(x, y)⟩ +
            φˢᶜₙₘ(x, y) * ⟨φˢᶜₙₘ(x, y), u₀(x, y)⟩ / ⟨φˢᶜₙₘ(x, y), φˢᶜₙₘ(x, y)⟩ +
            φˢˢₙₘ(x, y) * ⟨φˢˢₙₘ(x, y), u₀(x, y)⟩ / ⟨φˢˢₙₘ(x, y), φˢˢₙₘ(x, y)⟩
        ), where
    φᶜᶜₙₘ(x, y) = cos(2 * π * n * x / l) * cos(2 * π * m * y / l),
    φᶜˢₙₘ(x, y) = cos(2 * π * n * x / l) * sin(2 * π * m * y / l),
    φˢᶜₙₘ(x, y) = sin(2 * π * n * x / l) * cos(2 * π * m * y / l),
    φˢˢₙₘ(x, y) = sin(2 * π * n * x / l) * sin(2 * π * m * y / l), and
    λₙₘ = (2 * π / l)^2 * (n^2 + m^2) * α.
Note that the inner product of two functions g(x, y) and h(x, y) is defined as
    ⟨g(x, y), h(x, y)⟩ = ∫_0^l ∫_0^l g(x, y) h(x, y) dx dy.
When n = 0 or m = 0, some of the inner product denominators above are 0, but
this doesn't actually matter because the corresponding numerators are also 0 and
those terms can just be ignored.
So, the solution operator for the homogeneous PDE (with f = 0) is
    F(u₀)(x, y, t) =
        ∑_{n = 0}^∞ ∑_{m = 0}^∞ exp(-λₙₘ * t) * (
            φᶜᶜₙₘ(x, y) * ⟨φᶜᶜₙₘ(x, y), u₀(x, y)⟩ / ⟨φᶜᶜₙₘ(x, y), φᶜᶜₙₘ(x, y)⟩ +
            φᶜˢₙₘ(x, y) * ⟨φᶜˢₙₘ(x, y), u₀(x, y)⟩ / ⟨φᶜˢₙₘ(x, y), φᶜˢₙₘ(x, y)⟩ +
            φˢᶜₙₘ(x, y) * ⟨φˢᶜₙₘ(x, y), u₀(x, y)⟩ / ⟨φˢᶜₙₘ(x, y), φˢᶜₙₘ(x, y)⟩ +
            φˢˢₙₘ(x, y) * ⟨φˢˢₙₘ(x, y), u₀(x, y)⟩ / ⟨φˢˢₙₘ(x, y), φˢˢₙₘ(x, y)⟩
        ).
We can express the initial condition of our PDE using in terms of its Fourier
series as
    u₀(x, y) =
        ∑_{n = 0}^∞ ∑_{m = 0}^∞ (
            u₀ᶜᶜₙₘ * φᶜᶜₙₘ(x, y) +
            u₀ᶜˢₙₘ * φᶜˢₙₘ(x, y) +
            u₀ˢᶜₙₘ * φˢᶜₙₘ(x, y) +
            u₀ˢˢₙₘ * φˢˢₙₘ(x, y)
        ).
Now, consider the inhomogeneous PDE for which
    f(x, y, t) = f̂(t)(x, y) =
        ∑_{n = 0}^∞ ∑_{m = 0}^∞ exp(-λ′ₙₘ * t) * (
            f₀ᶜᶜₙₘ * φᶜᶜₙₘ(x, y) +
            f₀ᶜˢₙₘ * φᶜˢₙₘ(x, y) +
            f₀ˢᶜₙₘ * φˢᶜₙₘ(x, y) +
            f₀ˢˢₙₘ * φˢˢₙₘ(x, y)
        ).
(Note that, if we allow the time-dependence to have a form that is not
exponential, the resulting solution will contain non-elementary integrals.)
Duhamel's formula tells us that the solution to the inhomogeneous PDE is
    u(x, y, t) = F(u₀)(x, y, t) + ∫_0^t F(f̂(τ))(x, y, t - τ) dτ =
        ∑_{n = 0}^∞ ∑_{m = 0}^∞ exp(-λₙₘ * t) * (
            u₀ᶜᶜₙₘ * φᶜᶜₙₘ(x, y) +
            u₀ᶜˢₙₘ * φᶜˢₙₘ(x, y) +
            u₀ˢᶜₙₘ * φˢᶜₙₘ(x, y) +
            u₀ˢˢₙₘ * φˢˢₙₘ(x, y)
        ) +
        ∫_0^t ∑_{n = 0}^∞ ∑_{m = 0}^∞ (
            exp(-λ′ₙₘ * τ - λₙₘ * (t - τ)) * (
                f₀ᶜᶜₙₘ * φᶜᶜₙₘ(x, y) +
                f₀ᶜˢₙₘ * φᶜˢₙₘ(x, y) +
                f₀ˢᶜₙₘ * φˢᶜₙₘ(x, y) +
                f₀ˢˢₙₘ * φˢˢₙₘ(x, y)
            )
        ) dτ =
        ∑_{n = 0}^∞ ∑_{m = 0}^∞ exp(-λₙₘ * t) * (
            u₀ᶜᶜₙₘ * φᶜᶜₙₘ(x, y) +
            u₀ᶜˢₙₘ * φᶜˢₙₘ(x, y) +
            u₀ˢᶜₙₘ * φˢᶜₙₘ(x, y) +
            u₀ˢˢₙₘ * φˢˢₙₘ(x, y)
        ) +
        ∑_{n = 0}^∞ ∑_{m = 0}^∞ (
            exp(-λₙₘ * t) / (λₙₘ - λ′ₙₘ) * (exp((λₙₘ - λ′ₙₘ) * t) - 1)
        ) * (
            f₀ᶜᶜₙₘ * φᶜᶜₙₘ(x, y) +
            f₀ᶜˢₙₘ * φᶜˢₙₘ(x, y) +
            f₀ˢᶜₙₘ * φˢᶜₙₘ(x, y) +
            f₀ˢˢₙₘ * φˢˢₙₘ(x, y)
        ).
If we let λ′ₙₘ = λₙₘ + Δλₙₘ, this simplifies to
    u(x, y, t) =
        ∑_{n = 0}^∞ ∑_{m = 0}^∞ exp(-λₙₘ * t) * (
            u₀ᶜᶜₙₘ * φᶜᶜₙₘ(x, y) +
            u₀ᶜˢₙₘ * φᶜˢₙₘ(x, y) +
            u₀ˢᶜₙₘ * φˢᶜₙₘ(x, y) +
            u₀ˢˢₙₘ * φˢˢₙₘ(x, y)
        ) +
        ∑_{n = 0}^∞ ∑_{m = 0}^∞ exp(-λₙₘ * t) / Δλₙₘ * (1 - exp(-Δλₙₘ * t)) * (
            f₀ᶜᶜₙₘ * φᶜᶜₙₘ(x, y) +
            f₀ᶜˢₙₘ * φᶜˢₙₘ(x, y) +
            f₀ˢᶜₙₘ * φˢᶜₙₘ(x, y) +
            f₀ˢˢₙₘ * φˢˢₙₘ(x, y)
        ).
For the test case below, we will only use the φˢˢₙₘ eigenfunction for specific
values of n and m. In other words, we will pick some constants n, m, u₀, f₀, and
Δλ, and we will set
    u(x, y, 0) = u₀ * φˢˢₙₘ(x, y) and
    f(x, y, t) = f₀ * exp(-(λₙₘ + Δλ) * t) * φˢˢₙₘ(x, y).
We should then end up with the solution
    u(x, y, t) =
        (u₀ + f₀ / Δλ * (1 - exp(-Δλ * t))) * exp(-λₙₘ * t) * φˢˢₙₘ(x, y).
In addition, we will use nondimensionalization to replace our variables with
    x̂ = x / l,
    ŷ = y / l,
    t̂ = t / (l^2 / α),
    û(x̂, ŷ, t̂) = u(x, y, t) / u₀,
    f̂(x̂, ŷ, t̂) = f(x, y, t) / (u₀ / (l^2 / α)).
Note that this converts the time t into the "Fourier number" α * t / l^2.
We will then define the nondimensionalized constants
    λ̂ₙₘ = λₙₘ * l^2 / α = (2 * π)^2 * (n^2 + m^2),
    Δλ̂ = Δλ * l^2 / α, and
    f̂₀ = f₀ / (u₀ / (l^2 / α))
We will also rewrite the eigenfunction in terms of the new variables as
    φ̂ˢˢₙₘ(x̂, ŷ) = φˢˢₙₘ(x̂ * l, ŷ * l) = sin(2 * π * n * x̂) * sin(2 * π * m * ŷ).
Our simplified PDE then becomes
    ∂û/∂t̂ = Δû + f̂, where
    û(x̂, ŷ, 0) = φ̂ˢˢₙₘ(x̂, ŷ) and
    f̂(x̂, ŷ, t̂) = f̂₀ * exp(-(λ̂ₙₘ + Δλ̂) * t̂) * φ̂ˢˢₙₘ(x̂, ŷ).
Our solution then becomes
    û(x̂, ŷ, t̂) =
        (1 + f̂₀ / Δλ̂ * (1 - exp(-Δλ̂ * t̂))) * exp(-λ̂ₙₘ * t̂) * φ̂ˢˢₙₘ(x̂, ŷ).
In order to improve readability, we will drop the hats from all variable names.
```
