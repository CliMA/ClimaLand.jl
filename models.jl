struct LandRT
    a_PAR_w1::Any
    a_NIR_w1::Any
    a_PAR_b1::Any
    a_NIR_b1::Any
    a_PAR_w2::Any
    a_NIR_w2::Any
    a_PAR_b2::Any
    a_NIR_b2::Any
    sum_fvh::Any
    fh::Any
    Kdiff::Any
end
Flux.@layer LandRT
function reflected(α_leaf, α_g, T)
    α_c = α_leaf * (1 - T)# As LAI -> ∞, T -> 0, canopy albedo approaches leaf albedo. As LAI -> 0, α_c -> 0
    gs_sum = T/(1-α_c*α_g)
    R = α_c + T*α_g*gs_sum
    S_abs = (1-α_g)*gs_sum
    C_abs = (1-α_c-T)*(α_g*gs_sum+1) # How to ensure that this is positive?
    return R
end

function (m::LandRT)(x̂, LAI, μ, α_PAR_g, α_NIR_g, SW_d, f_diff)
    α_PAR = sigmoid(m.a_PAR_w2 * leakyrelu.(m.a_PAR_w1 * x̂ .+ m.a_PAR_b1) .+ m.a_PAR_b2)
    α_NIR = sigmoid(m.a_NIR_w2 * leakyrelu.(m.a_NIR_w1 * x̂ .+ m.a_NIR_b1) .+ m.a_NIR_b2)
    sum_fvh = sigmoid(m.sum_fvh)
    fh = sigmoid(m.fh) .* sum_fvh
    G = @. fh* μ + (1-sum_fvh) * 1/2 + (sum_fvh - fh) * sin(acos(μ))
    K_dir = min.(G ./ μ, 100.0f0)
    K_diff = relu(m.Kdiff)
    transmitted_fraction = @. exp(-K_dir * LAI) * (1 - f_diff) + f_diff * exp(- K_diff * LAI)
    upwelling_from_land_par = @. reflected(α_PAR, α_PAR_g, transmitted_fraction)
    upwelling_from_land_nir = @. reflected(α_NIR, α_NIR_g, transmitted_fraction)
    return @. (upwelling_from_land_nir + upwelling_from_land_par)/2
end


struct LandRT2Stream
    a_PAR_w1::Any
    a_NIR_w1::Any
    a_PAR_b1::Any
    a_NIR_b1::Any
    a_PAR_w2::Any
    a_NIR_w2::Any
    a_PAR_b2::Any
    a_NIR_b2::Any
    τ_PAR::Any
    τ_NIR::Any
    sum_fvh::Any
    fh::Any
end
function reflected_twostream(
    G::FT,
    α_leaf::FT,
    τ_leaf::FT,
    LAI::FT,
    cosθs::FT,
    α_soil::FT,
    frac_diff::FT,
) where {FT}
    Ω = 1.0f0
    α_soil = max(eps(FT), α_soil) # this prevents division by zero, below.
    cosθs = max(eps(FT), cosθs) # The insolations package returns cosθs=0 (nighttime), but this assumes cosθs >0
    K = min.(G ./ cosθs, 100.0f0)
    # Compute μ̄, the average inverse diffuse optical length per LAI
    μ̄ = 1 / max(2G, eps(FT))

    # Clip this to prevent dividing by zero; the sum must also be < 1. Note that using eps(FT)
    # as the threshold leads to numerical errors, so we use 1e-4
    ω = min(max(α_leaf + τ_leaf, FT(1e-4)), 1 - FT(1e-4))

    # Compute aₛ, the single scattering albedo
    aₛ =
        0.5 * ω * (1 - cosθs * log((abs(cosθs) + 1) / max(abs(cosθs), eps(FT))))

    # Compute β₀, the direct upscattering parameter
    β₀ = (1 / ω) * aₛ * (1 + μ̄ * K) / (μ̄ * K)

    # Compute β, the diffuse upscattering parameter
    diff = α_leaf - τ_leaf
    # With uniform distribution, Dickinson integral becomes following:
    c²θ̄ = pi * G / 4
    β = 0.5 * (ω + diff * c²θ̄) / ω

    # Compute coefficients for two-stream solution
    b = 1 - ω + ω * β
    c = ω * β
    d = ω * β₀ * μ̄ * K
    f = ω * μ̄ * K * (1 - β₀)
    h = √(b^2 - c^2) / μ̄
    σ = (μ̄ * K)^2 + c^2 - b^2

    u₁ = b - c / α_soil
    u₂ = b - c * α_soil
    u₃ = f + c * α_soil

    s₁ = exp(-h * LAI * Ω)
    s₂ = exp(-K * LAI * Ω)

    p₁ = b + μ̄ * h
    p₂ = b - μ̄ * h
    p₃ = b + μ̄ * K
    p₄ = b - μ̄ * K

    d₁ = p₁ * (u₁ - μ̄ * h) / s₁ - p₂ * (u₁ + μ̄ * h) * s₁
    d₂ = (u₂ + μ̄ * h) / s₁ - (u₂ - μ̄ * h) * s₁

    # h coefficients for direct upward flux
    h₁ = -d * p₄ - c * f
    h₂ =
        1 / d₁ * (
            (d - h₁ / σ * p₃) * (u₁ - μ̄ * h) / s₁ -
            p₂ * s₂ * (d - c - h₁ / σ * (u₁ + μ̄ * K))
        )
    h₃ =
        -1 / d₁ * (
            (d - h₁ / σ * p₃) * (u₁ + μ̄ * h) * s₁ -
            p₁ * s₂ * (d - c - h₁ / σ * (u₁ + μ̄ * K))
        )

    # h coefficients for direct downward flux
    h₄ = -f * p₃ - c * d
    h₅ =
        -1 / d₂ *
        (h₄ * (u₂ + μ̄ * h) / (σ * s₁) + (u₃ - h₄ / σ * (u₂ - μ̄ * K)) * s₂)
    h₆ =
        1 / d₂ *
        (h₄ / σ * (u₂ - μ̄ * h) * s₁ + (u₃ - h₄ / σ * (u₂ - μ̄ * K)) * s₂)

    # h coefficients for diffuse upward flux
    h₇ = c * (u₁ - μ̄ * h) / (d₁ * s₁)
    h₈ = -c * s₁ * (u₁ + μ̄ * h) / d₁

    # h coefficients for diffuse downward flux
    h₉ = (u₂ + μ̄ * h) / (d₂ * s₁)
    h₁₀ = -s₁ * (u₂ - μ̄ * h) / d₂
    n_layers = 20
    # Compute the LAI per layer for this canopy
    Lₗ = LAI / n_layers

    # Initialize the fraction absorbed value and layer counter
    F_abs = 0
    i = 0

    # Total light reflected from top of canopy
    F_refl = 0

    # Total light transmitted through the canopy on the downward pass
    F_trans = 0

    # Intialize vars to save computed fluxes from each layer for the next layer
    I_dir_up_prev = 0
    I_dir_dn_prev = 0
    I_dif_up_prev = 0
    I_dif_dn_prev = 0


    # Compute F_abs in each canopy layer
    while i <= n_layers

        # Compute cumulative LAI at this layer
        L = i * Lₗ

        # Compute the direct fluxes into/out of the layer
        I_dir_up =
            h₁ * exp(-K * L * Ω) / σ +
            h₂ * exp(-h * L * Ω) +
            h₃ * exp(h * L * Ω)
        I_dir_dn =
            h₄ * exp(-K * L * Ω) / σ +
            h₅ * exp(-h * L * Ω) +
            h₆ * exp(h * L * Ω)

        # Add collimated radiation to downward flux
        I_dir_dn += exp(-K * L * Ω)

        # Compute the diffuse fluxes into/out of the layer
        I_dif_up = h₇ * exp(-h * L * Ω) + h₈ * exp(h * L * Ω)
        I_dif_dn = h₉ * exp(-h * L * Ω) + h₁₀ * exp(h * L * Ω)

        # Energy balance giving radiation absorbed in the layer
        if i == 0
            I_dir_abs = 0
            I_dif_abs = 0
        else
            I_dir_abs = I_dir_up - I_dir_up_prev - I_dir_dn + I_dir_dn_prev
            I_dif_abs = I_dif_up - I_dif_up_prev - I_dif_dn + I_dif_dn_prev
        end

        if i == 1 # prev = top of layer
            F_refl =
                (1 - frac_diff) * I_dir_up_prev + (frac_diff) * I_dif_up_prev
        end
        if i == n_layers # not prev = bottom of layer
            F_trans = (1 - frac_diff) * I_dir_dn + (frac_diff) * I_dif_dn
        end


        # Add radiation absorbed in the layer to total absorbed radiation
        F_abs += (1 - frac_diff) * I_dir_abs + (frac_diff) * I_dif_abs

        # Save input/output values to compute energy balance of next layer
        I_dir_up_prev = I_dir_up
        I_dir_dn_prev = I_dir_dn
        I_dif_up_prev = I_dif_up
        I_dif_dn_prev = I_dif_dn

        # Move on to the next layer
        i += 1
    end

    # Reflected is reflected from the land surface. F_abs
    # refers to radiation absorbed by the canopy on either pass.
    # F_trans refers to the downward pass only. Note that
    # (1-α_soil)*F_trans + F_abs = total absorbed fraction by land, so the following
    # must hold
    # @assert (1 - α_soil) * FT(F_trans) + FT(F_abs) + FT(F_refl) ≈ 1
    # This is tested in test/standalo
    return FT(F_refl)
end
Flux.@layer LandRT2Stream 
function (m::LandRT2Stream)(x̂, LAI, μ, α_PAR_g, α_NIR_g, SW_d, f_diff)
    α_PAR = sigmoid(m.a_PAR_w2 * leakyrelu.(m.a_PAR_w1 * x̂ .+ m.a_PAR_b1) .+ m.a_PAR_b2)
    α_NIR = sigmoid(m.a_NIR_w2 * leakyrelu.(m.a_NIR_w1 * x̂ .+ m.a_NIR_b1) .+ m.a_NIR_b2)
    sum_fvh = sigmoid(m.sum_fvh)
    τ_NIR = @. sigmoid(m.τ_NIR) * (1- α_NIR)
    τ_PAR = @. sigmoid(m.τ_PAR)* (1- α_PAR)
    fh = sigmoid(m.fh) .* sum_fvh
    G = @. fh* μ + (1-sum_fvh) * 1/2 + (sum_fvh - fh) * sin(acos(μ))
    upwelling_from_land_par = @. reflected_twostream(G, α_PAR, τ_PAR, LAI, μ, α_PAR_g, f_diff)
    upwelling_from_land_nir = @. reflected_twostream(G, α_NIR, τ_NIR, LAI, μ, α_NIR_g, f_diff)
    return @. (upwelling_from_land_nir + upwelling_from_land_par)/2
end


function compute_r2(y, y_hat)
    ss_res = sum((y .- y_hat) .^ 2)
    ss_tot = sum((y .- mean(y)) .^ 2)
    return 1.0 - (ss_res / ss_tot)
end

function get_α_PAR(m, x̂)
    return sigmoid(m.a_PAR_w2 * leakyrelu.(m.a_PAR_w1 * x̂ .+ m.a_PAR_b1) .+ m.a_PAR_b2)
end

function get_α_NIR(m, x̂)
    return sigmoid(m.a_NIR_w2 * leakyrelu.(m.a_NIR_w1 * x̂ .+ m.a_NIR_b1) .+ m.a_NIR_b2)
end
