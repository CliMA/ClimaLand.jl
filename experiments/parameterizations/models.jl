struct SoilOnlyAlbedo{A <: Chain, B <: Chain, C <: Chain}
    wet_0::A
    dry_0::B
    logΔα::C
end
function SoilOnlyAlbedo(nfeatures_soil)
    dry_0 = Chain(
        Dense(nfeatures_soil => 1, sigmoid),
    )
    wet_0 = Chain(
        Dense(nfeatures_soil => 1, sigmoid),
    )
    logΔα = Chain(Dense(1 => 1),)
    args = (wet_0, dry_0, logΔα)
    return SoilOnlyAlbedo{typeof.(args)...}(args...)
end
Flux.@layer SoilOnlyAlbedo
Flux.trainable(model::SoilOnlyAlbedo) = (; model.dry_0, model.wet_0, model.logΔα)
function (model::SoilOnlyAlbedo)(x; soil_ids, μ_id, θ_id)
    μ = x[μ_id,:]
    θ = x[θ_id, :]
    dry_albedo = model.dry_0(x[soil_ids,:])
    wet_albedo = model.wet_0(x[soil_ids, :])
    α_soil = @. dry_albedo * (1 - θ) + wet_albedo * θ
    α_soil_zenith_corrected = min.(max.(α_soil .+ exp.(model.logΔα(μ)), Float32(0)), Float32(1))# α = α_0 + Δαe^(kμ)
    return α_soil_zenith_corrected
end


struct TotalAlbedoSingleSoilCLMCanopy{A <: Chain, B <: Chain, C <: Chain}
    α_soil_PAR::A
    α_soil_NIR::B
    logΔα::C
end

function TotalAlbedoSingleSoilCLMCanopy(nfeatures)
    α_soil_PAR = Chain(
        Dense(nfeatures => 1, sigmoid),
    )
    α_soil_NIR = Chain(
        Dense(nfeatures => 1, sigmoid),
    )
    logΔα = Chain(Dense(1 => 1),)
    args = (α_soil_PAR, α_soil_NIR, logΔα)
    return TotalAlbedoSingleSoilCLMCanopy{typeof.(args)...}(args...)
end
Flux.@layer TotalAlbedoSingleSoilCLMCanopy
Flux.trainable(model::TotalAlbedoSingleSoilCLMCanopy) = (; model.α_soil_PAR, model.α_soil_NIR, model.logΔα,)
function (model::TotalAlbedoSingleSoilCLMCanopy)(x; soil_ids, μ_id, frac_diff_id, LAI_id, α_PAR_id, α_NIR_id, τ_PAR_id, τ_NIR_id)
    μ = x[μ_id,:]
    frac_diff = x[frac_diff_id, :]
    LAI = x[LAI_id, :]
    α_soil = model.α_soil_PAR(x[soil_ids,:])
    α_soil_zenith_corrected = min.(max.(α_soil .+ exp.(model.logΔα(μ)), Float32(0)), Float32(1))# α = α_0 + Δαe^(kμ)
    τ_leaf = x[τ_PAR_id,:]
    α_leaf = x[α_PAR_id,:]
    yhat_PAR = simpler_canopy_sw_rt_two_stream.(α_leaf, τ_leaf, LAI, μ, α_soil_zenith_corrected, frac_diff) # reflected fraction, i.e. albedo

    α_soil = model.α_soil_NIR(x[soil_ids,:])
    α_soil_zenith_corrected = min.(max.(α_soil .+ exp.(model.logΔα(μ)), Float32(0)), Float32(1))# α = α_0 + Δαe^(kμ)
    τ_leaf = x[τ_NIR_id,:]
    α_leaf = x[α_NIR_id,:]
    yhat_NIR = simpler_canopy_sw_rt_two_stream.(α_leaf, τ_leaf, LAI, μ, α_soil_zenith_corrected, frac_diff) # reflected fraction, i.e. albedo
    return (yhat_PAR .+ yhat_NIR) ./ 2
end





struct TotalAlbedo{A <: Chain, B <: Chain, C <: Chain, D <: Chain,  E <: Chain}
    wet_0::A
    dry_0::B
    logΔα::C
    α_leaf::D
    sum_α_τ_leaf::E
end

function TotalAlbedo(nfeatures_soil, nfeatures_canopy)
    dry_0 = Chain(
        Dense(nfeatures_soil => 1, sigmoid),
    )
    wet_0 = Chain(
        Dense(nfeatures_soil => 1, sigmoid),
    )
    logΔα = Chain(Dense(1 => 1),)
    α_leaf = Chain(
        Dense(nfeatures_canopy => 1, sigmoid),
    )
    sum_α_τ_leaf = Chain(
        Dense(nfeatures_canopy => 1, sigmoid),
    )
    args = (wet_0, dry_0, logΔα, α_leaf, sum_α_τ_leaf)
    return TotalAlbedo{typeof.(args)...}(args...)
end
Flux.@layer TotalAlbedo
Flux.trainable(model::TotalAlbedo) = (; model.dry_0, model.wet_0, model.logΔα,model.α_leaf, model.sum_α_τ_leaf)
function (model::TotalAlbedo)(x; soil_ids, canopy_ids, μ_id, θ_id, frac_diff_id, LAI_id)
    μ = x[μ_id,:]
    frac_diff = x[frac_diff_id, :]
    LAI = x[LAI_id, :]
    θ = x[θ_id, :]
    dry_albedo = model.dry_0(x[soil_ids,:])
    wet_albedo = model.wet_0(x[soil_ids, :])
    α_soil = @. dry_albedo * (1 - θ) + wet_albedo * θ
    α_soil_zenith_corrected = min.(max.(α_soil .+ exp.(model.logΔα(μ)), Float32(0)), Float32(1))# α = α_0 + Δαe^(kμ)
    α_leaf = model.α_leaf(x[canopy_ids, :]) 
    sum_α_τ_leaf = model.sum_α_τ_leaf(x[canopy_ids, :])
    τ_leaf  = max.(sum_α_τ_leaf .- α_leaf,0.0f0)
    yhat = simpler_canopy_sw_rt_two_stream.(α_leaf, τ_leaf, LAI, μ, α_soil_zenith_corrected, frac_diff) # reflected fraction, i.e. albedo
    return yhat
end


struct TotalAlbedoSingleSoil{A <: Chain, C <: Chain, D <: Chain,  E <: Chain}
    α_soil::A
    logΔα::C
    ratio_α_τ_leaf::D
    sum_α_τ_leaf::E
end

function TotalAlbedoSingleSoil(nfeatures_soil, nfeatures_canopy)
    α_soil = Chain(
        Dense(nfeatures_soil => 1, sigmoid),
    )
    logΔα = Chain(Dense(1 => 1),)
    ratio_α_τ_leaf = Chain(
        Dense(nfeatures_canopy => 1, relu),
    )
    sum_α_τ_leaf = Chain(
        Dense(nfeatures_canopy => 1, sigmoid),
    )
    args = (α_soil, logΔα, ratio_α_τ_leaf, sum_α_τ_leaf)
    return TotalAlbedoSingleSoil{typeof.(args)...}(args...)
end
Flux.@layer TotalAlbedoSingleSoil
Flux.trainable(model::TotalAlbedoSingleSoil) = (; model.α_soil, model.logΔα,model.ratio_α_τ_leaf, model.sum_α_τ_leaf)
function (model::TotalAlbedoSingleSoil)(x; soil_ids, canopy_ids, μ_id, frac_diff_id, LAI_id)
    μ = x[μ_id,:]
    frac_diff = x[frac_diff_id, :]
    LAI = x[LAI_id, :]
    α_soil = model.α_soil(x[soil_ids,:])
    α_soil_zenith_corrected = min.(max.(α_soil .+ exp.(model.logΔα(μ)), Float32(0)), Float32(1))# α = α_0 + Δαe^(kμ)
    ratio_α_τ_leaf = model.ratio_α_τ_leaf(x[canopy_ids, :])
    sum_α_τ_leaf = model.sum_α_τ_leaf(x[canopy_ids, :])
    τ_leaf = sum_α_τ_leaf ./ (1 .+  ratio_α_τ_leaf)
    α_leaf = sum_α_τ_leaf .- τ_leaf
    yhat = simpler_canopy_sw_rt_two_stream.(α_leaf, τ_leaf, LAI, μ, α_soil_zenith_corrected, frac_diff) # reflected fraction, i.e. albedo
    return yhat
end

struct TotalAlbedoSingleSoilNoTau{A <: Chain,  B <: Chain, C <: Chain, D <: Chain}
    α_soil::A
    logΔα::B
    α_leaf::C
    LAI_sigmoid::D
end

function TotalAlbedoSingleSoilNoTau(nfeatures_soil, nfeatures_canopy)
    α_soil = Chain(
        Dense(nfeatures_soil => 1, sigmoid),
    )
    logΔα = Chain(Dense(1 => 1),)
    α_leaf = Chain(
        Dense(nfeatures_canopy => 1, sigmoid),
    )
    LAI_sigmoid = Chain(Dense(1 => 1), sigmoid)
    args = (α_soil, logΔα, α_leaf, LAI_sigmoid)
    return TotalAlbedoSingleSoilNoTau{typeof.(args)...}(args...)
end
Flux.@layer TotalAlbedoSingleSoilNoTau
Flux.trainable(model::TotalAlbedoSingleSoilNoTau) = (; model.α_soil, model.logΔα,model.α_leaf, model.LAI_sigmoid)
function (model::TotalAlbedoSingleSoilNoTau)(x; soil_ids, canopy_ids, μ_id, frac_diff_id, LAI_id)
    μ = x[μ_id,:]
    LAI = x[LAI_id, :]
    α_soil = model.α_soil(x[soil_ids,:])
    α_soil_zenith_corrected = min.(max.(α_soil .+ exp.(model.logΔα(μ)), Float32(0)), Float32(1))# α = α_0 + Δαe^(kμ)
    α_leaf = model.α_leaf(x[canopy_ids, :])
    α_leaf_zenith_corrected = min.(max.(α_leaf .+ exp.(model.logΔα(μ)), Float32(0)), Float32(1))# α = α_0 + Δαe^(kμ)

    sigmoid_lai = model.LAI_sigmoid(LAI)
    yhat = @. α_soil_zenith_corrected*(1-sigmoid_lai) + sigmoid_lai*α_leaf_zenith_corrected
    return yhat
end



function simpler_canopy_sw_rt_two_stream(
    α_leaf::FT,
    τ_leaf::FT,
    LAI::FT,
    cosθs::FT,
    α_soil::FT,
    frac_diff::FT;
    G = FT(0.5),
    Ω = FT(1),
    n_layers = UInt64(20)
) where {FT}
    
    K = G/max(cosθs, eps(Float32))
    
    
    # Compute μ̄, the average inverse diffuse optical length per LAI
    μ̄ = 1 / (2G)

    # Clip this to eps(Float32) to prevent dividing by zero
    ω = max(α_leaf + τ_leaf, eps(Float32))

    # Compute aₛ, the single scattering albedo
    aₛ = 0.5 * ω * (1 - cosθs * log((abs(cosθs) + 1) / abs(cosθs)))

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

    # Compute the LAI per layer for this canopy
    Lₗ = LAI / n_layers

    # Initialize the fraction absorbed value and layer counter
    F_abs = 0
    i = 0

    # Total light reflected form top of canopy
    F_refl = 0

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

        # Add collimated radiation to downard flux
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

        if i == 1
            F_refl = (1 - frac_diff) * I_dir_up + (frac_diff) * I_dif_up
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
    return FT(F_refl)
end


function put_on_grid(lat, lon, x)
    lat_bins = 150
    lon_bins = 360
    x_grid = zeros(lon_bins, lat_bins) .+ NaN
    for i in 1:lat_bins
        for j in 1:lon_bins
            mask = (lat .> -60 +(i-1)) .&& (lat .< -60 + i) .&& (lon .> -180 +(j-1)) .&& (lon .< -180 + j)
            if sum(mask) > 0
                x_grid[j,i] = mean(x[mask])
            end
        end
    end
    return Array(-59.5:1:89.5), Array(-179.5:1:179.5), x_grid
end


function put_on_grid_and_count(x, y)
    σx = std(x); μx = mean(x)
    σy = std(y); μy = mean(y)
    xmin = μx-σx*3; xmax = μx+σx*3; Δx = σx/4
    ymin = μy-σy*3; ymax = μy+σy*3; Δy = σy/4
    x_bins = xmin:Δx:xmax
    y_bins = ymin:Δy:ymax
    nbins = length(Array(x_bins))-1
    counts = zeros(nbins, nbins) .+ NaN
    for i in 1:nbins
        for j in 1:nbins
            mask = (x .> xmin +(i-1)*Δx) .&& (x .< xmin + i*Δx) .&& (y .> ymin +(j-1)*Δy) .&& (y .< ymin + j*Δy)
            counts[i,j] = sum(mask) ./ length(mask)
        end
    end
    return Array((x_bins[2:end] .+ x_bins[1:end-1]) ./2),Array((y_bins[2:end] .+ y_bins[1:end-1]) ./2), counts
end
