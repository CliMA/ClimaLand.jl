using CubedSphere
using CubedSphere.SphericalGeometry

using LinearAlgebra
using Statistics
using Random

"""
    geometric_spacing(N, ratio_raised_to_N_minus_one)

Construct a symmetric set of `N` face locations on the interval `[-1, 1]` with geometrically graded spacing away from
the domain center. Let `r` be the geometric ratio recovered from `ratio_raised_to_N_minus_one` via
`r = (ratio_raised_to_N_minus_one)^(1/(N-1))`. The gaps between consecutive faces grow by a factor of `r` as one moves
outward from the center, and the layout is mirrored about zero. Endpoints are fixed at `-1` and `+1`.

- **Odd `N`**: A face lies at `0`. Outward gaps are `Δx, Δx*r, Δx*r^2, …`, mirrored about `0`, and chosen so the
  rightmost face lands at `+1`.
- **Even `N`**: No face at `0`. Central faces are at `±Δx/2`; outward gaps are `Δx*r, Δx*r^2, …`, mirrored, and chosen
  so the rightmost face lands at `+1`.

# Arguments
- `N`: Number of faces (`N ≥ 2`).
- `ratio_raised_to_N_minus_one`: The value `r^(N-1)` for some geometric ratio `r > 0`. Values near `1` yield nearly
  uniform spacing. (Exactly `1` is not supported by the closed-form formulas used here.)

# Returns
- `faces`: A length-`N` monotonically increasing vector of face coordinates on `[-1, 1]` with geometric grading and
           symmetry: `faces[1] = -1`, `faces[N] = 1`, and `faces[i] = -faces[N+1-i]`.
"""
function geometric_spacing(N, ratio_raised_to_N_minus_one)
    ratio = ratio_raised_to_N_minus_one^(1/(N - 1))
    faces = zeros(N)

    if isodd(N)
        M = round(Int, (N + 1)/2)

        Δx = 1 * (ratio - 1) / (ratio^(M - 1) - 1)

        faces[M] = 0

        k = 0

        for i in M+1:N
            faces[i] = faces[i-1] + Δx * ratio^k
            faces[N+1-i] = -faces[i]
            k += 1
        end
    else
        M = Int(N/2)

        Δx = 1/((ratio^M - 1)/(ratio - 1) - 0.5)

        faces[M] = -0.5Δx
        faces[M+1] = 0.5Δx

        k = 1

        for i in M+2:N
            faces[i] = faces[i-1] + Δx * ratio^k
            faces[N+1-i] = -faces[i]
            k += 1
        end
    end

    faces[1] = -1
    faces[N] = 1

    return faces
end

"""
    exponential_spacing(N, k₀ByN)

Construct a symmetric set of `N` face locations on `[-1, 1]` with **exponentially graded** spacing away from the domain
center. Let `k₀ = k₀ByN * N`, and define an exponential map on the right half, `x(t) = α·exp(t/k₀) + β`, anchored so
that it passes through `(t₀, 0)` and `(t₁, 1)`, then mirror about zero. Endpoints are fixed at `-1` and `+1`.

- **Odd `N`** (`M = (N+1)/2`): a face lies at `0` (`faces[M] = 0`). The right-half faces use `t = 1, …, M` with
  `x(1) = 0`, `x(M) = 1`, and the left half is the negative mirror.
- **Even `N`** (`M = N/2`): no face at `0`. The two central faces straddle zero, with `0` midway between them; the
  right-half faces use `t = 2, …, M+1` anchored by `x(1.5) = 0`, `x(M+1) = 1`, and the left half is mirrored.

# Arguments
- `N`: Number of faces (`N ≥ 2`).
- `k₀ByN`: Grading parameter scaled by `N` (the code uses `k₀ = k₀ByN * N`). Larger `k₀ByN` yields spacing closer to
  uniform; smaller `k₀ByN` increases clustering near the center. Requires `k₀ByN > 0`.

# Returns
- `faces`: A length-`N` strictly increasing vector of face coordinates on `[-1, 1]` with symmetry
           `faces[i] = - faces[N+1-i]`, and end-points `faces[1] = -1`, `faces[N] = 1`.
"""
function exponential_spacing(N, k₀ByN)
    k₀ = k₀ByN * N
    faces = zeros(N)

    if isodd(N)
        M = round(Int, (N+1)/2)

        A = [exp(1/k₀) 1
             exp(M/k₀) 1]

        b = [0, 1]

        coefficients = A \ b

        faces[M:N] = coefficients[1] * exp.((1:M)/k₀) .+ coefficients[2]

        for i in 1:M-1
            faces[i] = -faces[N+1-i]
        end

        faces[M] = 0
    else
        M = Int(N/2)

        A = [exp(1.5/k₀)   1
             exp((M+1)/k₀) 1]

        b = [0, 1]

        coefficients = A \ b

        faces[M+1:N] = coefficients[1] * exp.((2:M+1)/k₀) .+ coefficients[2]

        for i in 1:M
            faces[i] = -faces[N+1-i]
        end
    end

    faces[1] = -1
    faces[N] = 1

    return faces
end

"""
    conformal_cubed_sphere_coordinates(Nx, Ny;
                                       spacing = UniformSpacing(),
                                       spacing_parameters = (; ratio_raised_to_Nx_minus_one = 10.5,
                                                               k₀ByNx = 0.45))

Generate computational-space coordinates `x` and `y` on `[-1, 1] × [-1, 1]` and map them
to Cartesian coordinates `(X, Y, Z)` on the sphere using `conformal_cubed_sphere_mapping`.
The arrays `X`, `Y`, and `Z` are of size `(Nx, Ny)` and correspond to `Nx × Ny` grid vertices,
defining `(Nx-1) × (Ny-1)` spherical quadrilateral cells of a conformal cubed sphere panel.

If `spacing = UniformSpacing()` then `x` and `y` have uniform spacing (equal increments),
so their tensor product defines a uniform grid on `[-1, 1] × [-1, 1]`. Otherwise,
`x` and `y` both have a symmetric graded spacing. Specifically:
- `spacing = GeometricSpacing()`: uses `geometric_spacing(N, ratio_raised_to_Nx_minus_one)` for each axis.
- `spacing = ExponentialSpacing()`: uses `exponential_spacing(N, k₀ByNx)` for each axis.

# Arguments
- `Nx, Ny`: Number of grid vertices along the `x` and `y` directions (≥ 2).

# Keyword arguments
- `spacing`: Either `UniformSpacing()` (default), `GeometricSpacing()`, or `ExponentialSpacing()`.
- `spacing_parameters`: Parameters required for various spacings.
                        Default: `(; ratio_raised_to_Nx_minus_one = 10.5, k₀ByNx = 0.45)`.
                        For `GeometricSpacing()`, `ratio_raised_to_Nx_minus_one` is interpreted
                        as `r^(Nx-1)`; for `ExponentialSpacing()`, `k₀ByNx` is used
                        as `k₀ = k₀ByNx * Nx`.

# Returns
- `x`, `y`: Computational-space vertex coordinates of lengths `Nx` and `Ny`.
- `X`, `Y`, `Z`: `(Nx, Ny)` arrays of Cartesian coordinates on the sphere, with
                 `X[i, j], Y[i, j], Z[i, j] = conformal_cubed_sphere_mapping(x[i], y[j])`.
"""
function conformal_cubed_sphere_coordinates(Nx, Ny;
                                            spacing = UniformSpacing(),
                                            spacing_parameters = (; ratio_raised_to_Nx_minus_one = 10.5,
                                                                    k₀ByNx = 0.45))

    x, y = cube_face_coordinates(spacing, Nx, Ny, params=spacing_parameters)

    X, Y, Z = cube_to_sphere(x, y)

    return x, y, X, Y, Z
end

function cube_face_coordinates(::UniformSpacing, Nx, Ny; params=nothing)
    x = range(-1, 1, length = Nx)
    y = range(-1, 1, length = Ny)

    return x, y
end

function cube_face_coordinates(::GeometricSpacing, Nx, Ny; params)
    # For Nx = Ny = 32 + 1, setting ratio = 1.0775 increases the minimum cell width by a factor of 1.92.
    # For Nx = Ny = 1024 + 1, setting ratio = 1.0042 increases the minimum cell width by a factor of 3.25.
    x = geometric_spacing(Nx, params.ratio_raised_to_Nx_minus_one)
    y = geometric_spacing(Ny, params.ratio_raised_to_Nx_minus_one)
    return x, y
end

function cube_face_coordinates(::ExponentialSpacing, Nx, Ny; params)
    # For Nx = Ny = 32 + 1, setting k₀ByNx = 15 increases the minimum cell width by a factor of 1.84.
    # For Nx = Ny = 1024 + 1, setting k₀ByNx = 10 increases the minimum cell width by a factor of 2.58.
    x = exponential_spacing(Nx, params.k₀ByNx)
    y = exponential_spacing(Ny, params.k₀ByNx)
    return x, y
end

"""
    specify_parameters(spacing)

Return a one-element vector `[θ]` containing the initial guess for the spacing parameters
used to build non-uniform conformal cubed sphere panels for different `spacing`s.

- `GeometricSpacing()`: uses a single parameter interpreted downstream as `ratio^(N-1)`
   for geometric grading.
- `ExponentialSpacing()`: uses a single parameter interpreted downstream as `k₀/N`
  for exponential grading.
"""
specify_parameters(::UniformSpacing) = nothing
specify_parameters(::GeometricSpacing) = [1.0775]
specify_parameters(::ExponentialSpacing) = [15]

"""
    specify_parameter_limits(spacing)

Return the lower and upper bounds for the single spacing parameter as a **2-element vector** `[min, max]`.
For a consistent interface with possible multi-parameter extensions, this vector is returned
wrapped in a one-element array: `[[min, max]]`. The bounds depend on `spacing`:

- `GeometricSpacing()` returns `[[min, max]]` for `ratio^(N-1)`.
- `ExponentialSpacing()` returns `[[min, max]]` for `k₀/N`.

Returns `[θ_limits]` where `θ_limits == [min, max]` for the single parameter in this implementation.
"""
specify_parameter_limits(::GeometricSpacing) = [[5, 15]]
specify_parameter_limits(::ExponentialSpacing) = [[0.4, 0.5]]

"""
    specify_random_parameters(N_ensemble, spacing)

Draw an ensemble of random parameter vectors within the limits from
`specify_parameter_limits(spacing)`. Each ensemble member is sampled
uniformly within its parameter bounds.

# Arguments
- `N_ensemble`: Number of ensemble members to generate.
- `spacing`: `GeometricSpacing()` or `ExponentialSpacing()`.

# Returns
- A vector of length `N_ensemble`, where each element is a one-element parameter vector `[θ]` (a single parameter in this
  implementation).
"""
function specify_random_parameters(N_ensemble, spacing)
    θ = specify_parameters(spacing)
    θ_limits = specify_parameter_limits(spacing)

    θᵣ = [[θ_limits[j][1] + (θ_limits[j][2] - θ_limits[j][1]) * rand() for j in 1:lastindex(θ)] for i in 1:N_ensemble]

    return θᵣ
end

"""
    specify_weights_for_model_diagnostics()

Return the weights applied to the model diagnostics used by the objective function.
The two diagnostics are, in order: `(1) normalized minimum cell width`,
`(2) deviation from isotropy`.

# Returns
- A 2-element vector of weights, e.g., `[10, 1]`.
"""
function specify_weights_for_model_diagnostics()
    weights = [10, 1]
    return weights
end

"""
    compute_model_diagnostics(X, Y, Z, minimum_reference_cell_area)

Compute the two model diagnostics for a given mapped grid:

1. **Normalized minimum cell width**: `sqrt(min(cell_area) / minimum_reference_cell_area)`,
   using spherical cell areas from `compute_cell_areas(X, Y, Z)` and a uniform-grid reference area.
2. **Deviation from isotropy**: the heuristic anisotropy measure from `compute_deviation_from_isotropy(X, Y, Z)`.

# Arguments
- `X, Y, Z`: `(Nx, Ny)` Cartesian coordinate arrays on the sphere.
- `minimum_reference_cell_area`: The minimum cell area of a reference (typically uniform) grid.

# Returns
- A 2-element vector `[normalized_minimum_cell_width, deviation_from_isotropy]`.
"""
function compute_model_diagnostics(X, Y, Z, minimum_reference_cell_area)
    cell_areas = compute_cell_areas(X, Y, Z)
    normalized_minimum_cell_width = sqrt(minimum(cell_areas)/minimum_reference_cell_area)

    deviation_from_isotropy = compute_deviation_from_isotropy(X, Y, Z)

    model_diagnostics = vcat(normalized_minimum_cell_width, deviation_from_isotropy)

    return model_diagnostics
end

"""
    compute_weighted_model_diagnostics(model_diagnostics)

Apply weights from `specify_weights_for_model_diagnostics()` to the unweighted diagnostics.

# Arguments
- `model_diagnostics`: A 2-element vector `[normalized_minimum_cell_width, deviation_from_isotropy]`.

# Returns
- A 2-element vector with weights applied, in the same order as the inputs.
"""
function compute_weighted_model_diagnostics(model_diagnostics)
    normalized_minimum_cell_width = model_diagnostics[1]
    deviation_from_isotropy = model_diagnostics[2]

    weights = specify_weights_for_model_diagnostics()

    weighted_model_diagnostics = vcat(weights[1] * normalized_minimum_cell_width, weights[2] * deviation_from_isotropy)

    return weighted_model_diagnostics
end

"""
    forward_map(Nx, Ny, spacing, θ)

Evaluate the (weighted) model diagnostics for a non-uniform conformal cubed-sphere panel defined by
parameters `θ`. The steps are:
1. Clamp `θ` to parameter limits from `specify_parameter_limits(spacing)`.
2. Build a **reference** uniform grid via `conformal_cubed_sphere_coordinates(Nx, Ny; spacing=UniformSpacing())`
   and compute the `minimum_reference_cell_area`.
3. Build the **non-uniform** grid via `conformal_cubed_sphere_coordinates(Nx, Ny; spacing, ...)`,
   passing the parameters in `θ` according to `spacing` type.
4. Compute model diagnostics and then apply weights.

# Arguments
- `Nx, Ny`: Number of grid vertices along the panel coordinates.
- `spacing`: `GeometricSpacing()` (interprets `θ[1]` as `ratio^(Nx-1)`) or
             `ExponentialSpacing()` (interprets `θ[1]` as `k₀/Nx`).
- `θ`: Parameter vector (one element in this implementation).

# Returns
- A 2-element vector of **weighted** model diagnostics.
"""
function forward_map(Nx, Ny, spacing, θ)
    θ_limits = specify_parameter_limits(spacing)

    for i in 1:lastindex(θ)
        θ[i] = clamp(θ[i], θ_limits[i][1], θ_limits[i][2])
    end

    x_reference, y_reference, X_reference, Y_reference, Z_reference =
        conformal_cubed_sphere_coordinates(Nx, Ny; spacing = UniformSpacing())

    cell_areas = compute_cell_areas(X_reference, Y_reference, Z_reference)
    minimum_reference_cell_area = minimum(cell_areas)

    spacing_parameters = (; ratio_raised_to_Nx_minus_one = θ[1], k₀ByNx = θ[1])
    x, y, X, Y, Z = conformal_cubed_sphere_coordinates(Nx, Ny; spacing, spacing_parameters)

    model_diagnostics = compute_model_diagnostics(X, Y, Z, minimum_reference_cell_area)

    weighted_model_diagnostics = compute_weighted_model_diagnostics(model_diagnostics)

    return weighted_model_diagnostics
end

"""
    specify_ideal_weighted_model_diagnostics()

Return the target (ideal) values for the **weighted** model diagnostics used by the inversion.
By default, the target normalized minimum cell width is `4` and the target deviation from
isotropy is `0`, then weights are applied in the same order.

# Returns
- A 2-element vector of target **weighted** diagnostics.
"""
function specify_ideal_weighted_model_diagnostics()
    normalized_minimum_cell_width = 4

    deviation_from_isotropy = 0

    weights = specify_weights_for_model_diagnostics()

    ideal_weighted_model_diagnostics = vcat(weights[1] * normalized_minimum_cell_width,
                                            weights[2] * deviation_from_isotropy)

    return ideal_weighted_model_diagnostics
end

"""
    optimize!(Nx, Ny, spacing, θ; N_iterations=10, Δt=1)

Run an Ensemble Kalman Inversion (EKI) to tune the spacing parameter(s) for a non-uniform conformal
cubed sphere panel. An ensemble of parameter vectors `θ` is iteratively updated so that the
**weighted** model diagnostics produced by `forward_map` match the ideal targets from
`specify_ideal_weighted_model_diagnostics()`.

At each iteration:
1. **Forward evaluations (parallelized):** For each ensemble member `θ[n]`, compute the predicted diagnostics
   `G[n] = forward_map(Nx, Ny, spacing, θ[n])`.
2. **Ensemble statistics:** Form means `θ̄ = mean(θ)` and `G̅ = mean(G)`, then compute
   - Cross-covariance `Cᵘᵖ = cov(θ, G)` (shape: `nθ × ndata`), and
   - Data covariance `Cᵖᵖ = cov(G, G)` (shape: `ndata × ndata`),
   using the unbiased `(N_ensemble-1)` denominator.
3. **Perturbed observations:** For each member, build `y[n] = ideal + Δt * η[n]` where `η[n] ~ N(0, I)`.
   Here `Δt` sets the **perturbation magnitude** of the observations.
4. **Residuals:** `r[n] = y[n] - G[n]`.
5. **Implicit update (Kalman-like step):** Update each parameter vector via `θ[n] ← θ[n] + K * r[n]`, with
   `K = Cᵘᵖ * (Cᵖᵖ + I/Δt)⁻¹`, implemented by solving the linear system with a Cholesky factorization of
   `Cᵖᵖ + I/Δt`. The same `Δt` also acts as an **implicit damping/step-size control**:
   smaller `Δt` ⇒ stronger regularization and smaller updates;
   larger `Δt` ⇒ weaker regularization and larger, noisier updates.
6. **Monitoring:** Report `error = ‖mean(r)‖` and store a snapshot of the ensemble.

This function **mutates** the input ensemble `θ` in place (it becomes the final ensemble) and records the full ensemble
after each iteration.

# Arguments
- `Nx, Ny`: Number of grid vertices along the panel coordinates.
- `spacing`: `GeometricSpacing()` or `ExponentialSpacing()`.
- `θ`: Initial ensemble — a vector of one-element parameter vectors `[θ]` (length `N_ensemble`).

# Keyword Arguments
- `N_iterations`: Number of EKI iterations. Default: 10
- `Δt`: Pseudo-time step that sets both the observation perturbation scale and the implicit damping
        in the linear solve `(Cᵖᵖ + I/Δt)`.

# Returns
- `θ_series`: A length `N_iterations + 1` vector; each entry is a **snapshot** of the full ensemble
  (index 1 is the initial ensemble; the last is the final ensemble). The input `θ` is also mutated
  to the final state.
"""
function optimize!(Nx, Ny, spacing, θ;
                   N_iterations = 10,
                   Δt = 1,
                   verbose = false)
    ideal_data = specify_ideal_weighted_model_diagnostics()
    model_data = forward_map(Nx, Ny, spacing, mean(θ))

    nData = length(ideal_data)
    N_ensemble = length(θ)

    θ_series = [copy(θ)]

    error = norm(model_data - ideal_data)

    if verbose
        @info("\nIteration 0 with error $error")
    end

    G = [copy(model_data) for i in 1:N_ensemble]

    # EKI iteration is equivalent to a time step of the above equation.
    @inbounds for i in 1:N_iterations
        θ̄ = mean(θ)

        # Evaluating the forward map for all ensemble members. This is the most expensive step because it needs to run
        # the model N_ensemble times. For the moment our model is simple, but imagine doing this with a full climate
        # model! Luckily this step is embarassingly parallelizeable.
        Threads.@threads for n in 1:N_ensemble
            G[n] .= forward_map(Nx, Ny, spacing, θ[n]) # Error handling needs to go here.
        end

        # The ensemble mean output of the models
        G̅ = mean(G)

        # Calculating the covariances to be used in the update steps
        Cᵘᵖ = (θ[1] - θ̄) * (G[1] - G̅)'
        Cᵖᵖ = (G[1] - G̅) * (G[1] - G̅)'

        for j = 2:N_ensemble
            Cᵘᵖ += (θ[j] - θ̄) * (G[j] - G̅)'
            Cᵖᵖ += (G[j] - G̅) * (G[j] - G̅)'
        end

        Cᵘᵖ *= 1 / (N_ensemble - 1)
        Cᵖᵖ *= 1 / (N_ensemble - 1)

        # Ensemblize the data (adding the random noise η).
        y = [ideal_data + Δt * randn(nData) for i in 1:N_ensemble]

        # The residual from our observations
        r = y - G

        # Update the parameters using implicit pseudo-time-stepping, which involves solving a linear system.
        Cᵖᵖ_factorized = cholesky(Symmetric(Cᵖᵖ + 1 / Δt * LinearAlgebra.I))

        for j in 1:N_ensemble
            θ[j] .+= Cᵘᵖ * (Cᵖᵖ_factorized \ r[j])
        end

        error = norm(mean(r))
        if verbose
            @info "Iteration $i with error $error"
        end
        push!(θ_series, copy(θ))
    end

    return θ_series
end

"""
    optimized_non_uniform_conformal_cubed_sphere_coordinates(Nx, Ny, spacing)

High-level driver that uses EKI to optimize the non-uniform spacing parameter for a conformal
cubed sphere panel, then builds and returns the corresponding grid.

Procedure:
1. Create a random ensemble of parameters within limits (`N_ensemble = 40`, reproducible seed).
2. Run `optimize!` to fit the **weighted** diagnostics to their ideal targets.
3. Build the optimized grid with `conformal_cubed_sphere_coordinates(Nx, Ny; spacing, ...)` using the
   mean optimized parameter.

For `GeometricSpacing()`, the parameter is `ratio^(Nx-1)`; for `ExponentialSpacing()`, the parameter is `k₀/Nx`.

# Arguments
- `Nx, Ny`: Number of grid vertices along panel coordinates.
- `spacing`: `GeometricSpacing()` or `ExponentialSpacing()`.

# Returns
- `x, y`: Computational-space coordinates of lengths `Nx` and `Ny`.
- `X, Y, Z`: `(Nx, Ny)` Cartesian coordinates of the vertices of the **optimized** non-uniform
             conformal cubed sphere panel.
"""
function optimized_non_uniform_conformal_cubed_sphere_coordinates(Nx, Ny, spacing;
                                                                  verbose = false)
    N_ensemble = 40 # Choose N_ensemble to be at least 4 times the number of parameters.

    if verbose
        @info "Optimize non-uniform conformal cubed sphere for Nx = $Nx and Ny = $Ny"
    end

    begin
        Random.seed!(123)
        θᵣ = specify_random_parameters(N_ensemble, spacing)
        θᵢ = deepcopy(θᵣ)

        θ_series = optimize!(Nx, Ny, spacing, θᵣ; N_iterations = 10)
    end

    if spacing isa GeometricSpacing
        θ_name = "ratio_raised_to_Nx_minus_one"
    elseif spacing isa ExponentialSpacing
        θ_name = "k₀ByNx"
    end

    if verbose
        println("\nThe unoptimized parameters are: $θ_name = $(round(mean(θᵢ)[1], digits=2))\n")
        println("\nThe optimized parameters are: $θ_name = $(round(mean(θᵣ)[1], digits=2))\n")
    end

    spacing_parameters = (; ratio_raised_to_Nx_minus_one = mean(θᵣ)[1],
                            k₀ByNx = mean(θᵣ)[1])
    x, y, X, Y, Z =
        conformal_cubed_sphere_coordinates(Nx, Ny; spacing, spacing_parameters)

    return x, y, X, Y, Z
end

"""
    compute_deviation_from_isotropy(X, Y, Z; radius = 1)

Compute a scalar measure of the deviation from isotropy for a spherical grid (e.g., a conformal
cubed-sphere panel), defined by the Cartesian coordinate arrays `X`, `Y`, and `Z`. Each of
`X`, `Y`, and `Z` is a 2D array of size `(Nx, Ny)` holding the Cartesian coordinates of
the grid vertices on the sphere, such that the point at `(i, j)` corresponds to
`(X[i, j], Y[i, j], Z[i, j])`. The grid therefore contains `(Nx−1) × (Ny−1)` spherical
quadrilateral cells.

For each quadrilateral cell, the function computes the great-circle distances of its four edges
on the sphere, evaluates the sum of absolute differences between consecutive edge lengths as a
measure of cell anisotropy, and then returns the Euclidean norm of these deviations over the
entire grid.

# Arguments
- `X`, `Y`, `Z`: `(Nx, Ny)` arrays of Cartesian coordinates of grid vertices on the sphere.

# Keyword Arguments
- `radius`: Sphere radius. Default: 1.

# Returns
- A non-negative scalar quantifying the overall deviation from isotropy in the grid.
  Larger values correspond to more anisotropic grids.

# Examples

A regular latitude-longitude grid:

```jldoctest 1
using CubedSphere: compute_deviation_from_isotropy

Nx, Ny = 3, 3
lons = range(-π/4, π/4, length = Nx)
lats = range(-π/6, π/6, length = Ny)

X = [cos(φ) * cos(λ) for λ in lons, φ in lats]
Y = [cos(φ) * sin(λ) for λ in lons, φ in lats]
Z = [sin(φ)          for λ in lons, φ in lats]

compute_deviation_from_isotropy(X, Y, Z)

# output

1.6552138747243959
```

The vertices of a tetrahedron:

```jldoctest
using CubedSphere: compute_deviation_from_isotropy

V = [   0     0     1 ;
     2√2/3    0   -1/3;
     -√2/3  √6/3  -1/3;
     -√2/3 -√6/3  -1/3]

X = reshape(V[:, 1], 2, 2)
Y = reshape(V[:, 2], 2, 2)
Z = reshape(V[:, 3], 2, 2)

isapprox(compute_deviation_from_isotropy(X, Y, Z), 0, atol=1e-14)

# output

true
```
"""
function compute_deviation_from_isotropy(X, Y, Z; radius=1)
    Nx, Ny = size(X)
    deviation_from_isotropy = zeros(Nx-1, Ny-1)

    for j in 1:Ny-1, i in 1:Nx-1
        a₁, a₂, a₃, a₄ = spherical_quadrilateral_vertices(X, Y, Z, i, j)

        # Compute the arc lengths (distances) between the points a₁ and a₂,
        # a₂ and a₃, a₃ and a₄, and a₄ and a₁ on the sphere.
        d₁ = spherical_distance(a₁, a₂; radius)
        d₂ = spherical_distance(a₂, a₃; radius)
        d₃ = spherical_distance(a₃, a₄; radius)
        d₄ = spherical_distance(a₄, a₁; radius)

        # Compute the deviation from isotropy.
        deviation_from_isotropy[i, j] = abs(d₁ - d₂) + abs(d₂ - d₃) + abs(d₃ - d₄) + abs(d₄ - d₁)
    end

    return norm(deviation_from_isotropy)
end
