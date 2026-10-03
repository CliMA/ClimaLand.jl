# # Mixed implicit and explicit timestepping
# The goal of this tutorial is to describe how
# the timestepping of a ClimaLand model is carried out.
# We will use forward and backward Euler as a demonstration,
# but higher order methods are available in ClimaTimesteppers.

# # Explicit vs. implicit stepping
# Given a differential equation for a prognostic variable
# Y

# ```
# \frac{d Y }{dt} = g(Y, t) + h(Y, t),
# ```

# an explicit (forward) Euler step would entail

# ```
# Y(t+\Delta t) = Y(t) + [g(Y(t), t) + h(Y(t),t)] \times \Delta t,
# ```

# while an implicit (backward) Euler step would entail

# ```
# Y(t+\Delta t) = Y(t) + [g(Y(t+\Delta t), t) + h(Y(t+\Delta t),t)] \times \Delta t,
# ```

# which reqires us to solve an implicit equation for  Y(t+ Δt). We
# usually do so using Newton's method, which requires the derivative of
# the entire right hand side with respect to the variable we are solving
# for, `Y(t+ Δt)`. This is called the Jacobian.

# Sometimes certain terms must be stepped implicit for numerical stability,
# while others are more slowly varying or stable. In this case, a mixed
# approach would be

# ```
# Y(t+\Delta t) = Y(t) + [g(Y(t+\Delta t), t) + h(Y(t),t)] \times \Delta t,
# ```

# assuming that `h` is the slow term, and `g` is the fast term. Note that
# solving this implicit equation for `(Y(t+ Δt)` with Newton's method
# would be similar to that of the fully implicit approach, but with an
# approximated Jacobian (neglecting ∂h/∂Y).

# If our timestepping scheme involves evaluating all right-hand-side
# tendencies at the current (known) value of a prognostic variable,
# we refer to that prognostic variable as explicit. If any of the
# right-hand-side tendencies are evaluated at the next (unknown) value
# of the prognostic variable, we refer to it as implicit. In the latter
# case, the Jacobian would include a term like `∂ tendency/∂ variable`,
# even if it is an approximate (not exact) form of the Jacobian.

# # Implicit and explicit prognostic variables of the land model
# We treat two prognostic variables of the soil model (ϑ_l,
# ρe_int) and the canopy temperature implicitly, and 
# the canopy water content, the soil ice content,
# and all prognostic
# variables of the snow model explicitly.

# # Implicit vs explicit tendencies - not complete as of 4/23/25
# Implicit:
# - Vertical contribution of the divergence of the Darcy flux
# - Vertical contribution of the divergence of diffusive heat flux
# - Vertical contribution of the divergence of heat flux due to Darcy flow
# - Canopy temperature (except root extraction of energy)
# - SHF, LHF, evaporation, and sublimation of soil (note that these are explicit in θ_i!)
# - Soil radiation (does not contribute to Jacobian)
# - Subsurface runoff (this is computed in the same function the same time as surface runoff, but does not contribute to Jacobian.)

# Explicit
# - Phase changes in soil
# - Root extraction
# - all snow tendencies
# - Darcy flux within the canopy
