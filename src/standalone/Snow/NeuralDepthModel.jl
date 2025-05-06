using Flux
using ClimaLand
using ClimaLand.Snow:
    AbstractDensityModel, SnowModel, SnowParameters, snow_bulk_density
import ClimaLand.Snow:
    density_prog_vars,
    density_prog_types,
    density_prog_names,
    update_density_and_depth!,
    update_density_prog!
using Thermodynamics
using StaticArrays
import Adapt

export NeuralDepthModel

"""
    NeuralDepthModel{FT <: AbstractFloat} <: AbstractDensityModel{FT}
Establishes the density parameterization where snow density is calculated from a prognostic snow depth variable,
along with the prognostic SWE variable, using a neural network for the rate of change of snow depth, `dz/dt`.
The input to the network are temporally averaged, which we achieve using an exponentially moving average,
with a rate of `\alpha`.

During construction, all Arrays in `z_model` are converted to StaticArrays (this is done
for gpu compatibilty). For gpu compatibilty, `z_model` must not allocate or take the adjoint
of a [1xN] array, where N is the number of batched input vectors.
"""
struct NeuralDepthModel{FT, MD <: Flux.Chain} <: AbstractDensityModel{FT}
    "The Flux neural network for compute dz/dt"
    z_model::MD
    "The inverse of the averaging window time (1/s)"
    α::FT
    function NeuralDepthModel(z_model::MD, α::FT) where {FT, MD}
        static_z_model = Adapt.adapt(SArray, z_model)
        new{FT, typeof(static_z_model)}(static_z_model, α)
    end
end

"""
    get_network_weights(FT::DataType)
Retrieves the network weights for the snow depth model formulated in Charbonneau et. al., 2025,
and stores/returns them within a dictionary to enable easy loading of the weights into the generated
neural model.
Note that because we are loading from an unchanging document of pre-trained weights,
the weight sizes are hard-coded in the function below.

# Arguments
- `FT::DataType`: the type of float for which the the parameters are converted (Float64 or Float32)
"""
function get_network_weights(FT::DataType)
    weightfile =
        readlines(ClimaLand.Artifacts.neural_depth_model_weights_path())
    @assert length(weightfile) == 56 "Artifact path might need updating, expected 56 lines and got $(length(weightfile))"
    data = Dict()
    data["SCALING"] = parse.(FT, split(weightfile[2], " "))
    data["l1_WEIGHT"] =
        hcat([parse.(FT, split(weightfile[i], " ")) for i in 5:32]...)'
    data["l1_BIAS"] = parse.(FT, split(weightfile[35], " "))
    data["l2_WEIGHT"] =
        hcat([parse.(FT, split(weightfile[i], " ")) for i in 38:44]...)'
    data["l2_BIAS"] = parse.(FT, split(weightfile[47], " "))
    data["l3_WEIGHT"] = parse.(FT, split(weightfile[50], " "))'
    data["l3_BIAS"] = parse.(FT, split(weightfile[53], " "))
    data["FINALSCALE"] = parse(FT, strip(weightfile[56]))
    return data
end

"""
    connectfunc(pred::AbstractArray{FT}, input::AbstractArray{FT}) where {FT <: AbstractFloat}
Builds the connection function for the SkipConnection in the neural network formulated in
Charbonneau et. al., creating the upper and lower boundary conditions to apply to the predicted
output by the neural model.
Note that because we are loading a pre-trained model, the indices of certain variables are hard coded in the function below.
The first function refers to the upper boundary, which is max(0, pred) whenever snow precipitation is nonzero.
The lower boundary refers to ensuring nonnegativity of the snowpack height, which will get scaled by the model's timescale constant
in a separate layer.

# Arguments
- `pred::AbstractArray`: the array of predicted values by the network
- `input::AbstractArray`: the array of input features to the network to make the above predicted values.
"""
function connectfunc(
    pred::AbstractArray{FT},
    input::AbstractArray{FT},
)::AbstractArray{FT} where {FT <: AbstractFloat}
    return [
        @. (input[7, :]' > 0) * relu(pred)
        input[1, :]'
        pred
    ]
end

# second method to deal with adjoint of zero dimensional SArray
function connectfunc(
    pred::SVector{1, FT},
    input::SVector{7, FT},
)::AbstractArray{FT} where {FT <: AbstractFloat}
    return SVector((input[7] > 0) * relu(pred[1]), input[1], pred[1])
end


"""
    build_z_model(FT)
Create the neural network structure from the paper.
This does not set the predictive/scaling weight values - these must be loaded or the model can be
trained on additional data.
Note that because we are loading a pre-trained model, the number of features input, the size of the output, etc, are hard coded in the function below.

# Arguments
- `FT::DataType`: the type of float for which the the parameters are converted (Float64 or Float32)
"""
function build_z_model(FT::DataType)
    nfeatures = 7
    k = 4
    in_scales = FT.(ones(7))
    scale_and_relus = Matrix{FT}([1 0 0; 0 1 0; 1 0 -1])
    get_min = Matrix{FT}([1 1 -1; 0 1 0])
    final_mul = Matrix{FT}([1 -1])

    model = Chain(
        pred = SkipConnection(
            Chain(
                scale = x -> x .* in_scales,
                l1 = Dense(nfeatures, k * nfeatures, relu),
                l2 = Dense(k * nfeatures, nfeatures, elu),
                l3 = Dense(nfeatures, 1),
            ),
            connectfunc,
        ),
        apply_relus = Dense(scale_and_relus, false, relu),
        apply_upper = Dense(get_min, false, relu),
        sumlayer = Dense(final_mul, false),
    )
    if FT == Float32
        return Flux.f32(model)
    elseif FT == Float64
        return Flux.f64(model)
    else
        error(
            "Conversion of Flux Chain weights to the desired float-type is not supported. Please choose Float32 or Float64, or add the desired
            functionality to build_z_model() within src/standalone/snow/NeuralDepthModel.jl.",
        )
    end
end


"""
    setoutscale!(model, scale)
Set the physical scaling parameter for model usage (i.e. rectifying scaling done on model input).
This scaling constant is combined with a -1 factor that enforces the boundary for speed and memory purposes -
for this reason, we do not recommend first `extracting`/storing the weight this function refers to,
and then using that value with this function to `reset` the weight the intial/stored value, as this will have unintended
consequnces without accounting for the extra factor. Only use this function to scale the output
of the predictive component by the constant `scale`.

# Arguments
- `model::Chain`: the neural model to be used.
- `scale::Real`: the scaling parameter to return data to applicable units.
"""
function setoutscale!(model, scale::Real)
    FT = eltype(model[:apply_relus].weight)
    model[:apply_relus].weight[end, 3] = FT(-1 * scale)
end


"""
    settimescale!(model, dt)
Set the timescale parameter for model usage

# Arguments
- `model::Chain`: the neural model to be used.
- `dt::Real`: the number of seconds per timestep for usage.
"""
function settimescale!(model, dt::Real)
    FT = eltype(model[:apply_relus].weight)
    model[:apply_relus].weight[2, 2] = FT(1.0 / dt)
end

"""
    get_znetwork()
Return the snow-depth neural network from Charbonneau et al (2025; https://arxiv.org/abs/2412.06819),
and set the network's scaling such that snowpack depth remains nonnegative for a default timestep of 86400 seconds (1 day).
Note that because we are loading a pre-trained model, the number of features input, the size of the output, etc, are hard coded in the function below.

# Arguments
- `FT::DataType`: (Optional) the type of float for which the the parameters are converted (Float64 or Float32). Default is Float32
"""
function get_znetwork(; FT = Float32)
    network = build_z_model(FT)
    data = get_network_weights(FT)
    network[:pred].layers.layers[:scale].in_scales .= data["SCALING"]
    network[:pred].layers.layers[:l1].weight .= data["l1_WEIGHT"]
    network[:pred].layers.layers[:l1].bias .= data["l1_BIAS"]
    network[:pred].layers.layers[:l2].weight .= data["l2_WEIGHT"]
    network[:pred].layers.layers[:l2].bias .= data["l2_BIAS"]
    network[:pred].layers.layers[:l3].weight .= data["l3_WEIGHT"]
    network[:pred].layers.layers[:l3].bias .= data["l3_BIAS"]
    setoutscale!(network, data["FINALSCALE"])
    settimescale!(network, 86400)
    return network
end

"""
    NeuralDepthModel(FT::DataType; Δt::Union{Nothing, FT}=nothing; model::Flux.Chain, α::Union{FT, Nothing})
An outer constructor for the NeuralDepthModel density parameterization for usage in a snow model.

This model utilizes the neural network formulated in Charbonneau et. al. (2025) to predict the rate of change of snow height dz/dt.
Since the default network receives inputs representing the average value per feature per day (86400 s), an exponential-weighting
scheme is used to track the moving average of the input fields. A default weighting is supplied
as `\alpha = 2/86400 s⁻¹` following common practice of choosing a weighting of 2/(desired temporal window size),
though the user may specify their own value by keyword argument. Note that this default value is
unstable (whenever `1/\alpha > model timestep`) for model timesteps larger than half a day (43200 seconds).
The user should also specify a model timestep Δt in seconds which sets the model scaling in order to guarantee the snow depth remains nonnegative."
"""
function NeuralDepthModel(
    FT::DataType;
    Δt::Union{AbstractFloat, Nothing} = nothing,
    α::Union{AbstractFloat, Nothing} = nothing,
)
    usemodel = get_znetwork(FT = FT)
    weight = !isnothing(α) ? FT(α) : FT(2 / 86400)
    if !isnothing(Δt)
        settimescale!(usemodel, Δt) #assumes Δt is provided in seconds
        if (Δt > 43200) && isnothing(α)
            error("Please supply a weight for the
                   exponential moving average, the
                   default value is unstable in this case.")
        end
    end
    return NeuralDepthModel(usemodel, weight)
end

#Define the additional prognostic variables needed for using this parameterization:
ClimaLand.Snow.density_prog_vars(m::NeuralDepthModel) =
    (:Z, :P_avg, :T_avg, :R_avg, :Qrel_avg, :u_avg)
ClimaLand.Snow.density_prog_types(m::NeuralDepthModel{FT}) where {FT} =
    (FT, FT, FT, FT, FT, FT)
ClimaLand.Snow.density_prog_names(m::NeuralDepthModel) =
    (:surface, :surface, :surface, :surface, :surface, :surface)


#Extend/Define the appropriate functions needed for this parameterization:
"""
    update_density_and_depth!(ρ_snow, z_snow,density::NeuralDepthModel, Y, p, params::SnowParameters,)

Updates the snow density and depth in place given the current model state. Extends the update_density_and_depth!
function for the NeuralDepthModel type.
"""
function update_density_and_depth!(
    ρ_snow,
    z_snow,
    density::NeuralDepthModel,
    Y,
    p,
    params::SnowParameters,
)
    @. z_snow = Y.snow.Z
    @. ρ_snow = snow_bulk_density(Y.snow.S, z_snow, params)
end


"""
    swe_snow_area(S::FT, scf::FT, z::FT)::FT where {FT}
Estimates the SWE over the snow-covered portion of a grid cell, assming Z is provided
as the depth over the snow-covered portion and S is provided as the depth over the whole cell area.
Used in computing the time derivative of snow depth.
"""
function swe_snow_area(S::FT, scf::FT, z::FT)::FT where {FT}
    min_scf = 0.05 #placeholder until Katherine and Andy discuss how to implement this
    return min(z, S / max(scf, min_scf)) #double check if this is the formula Katherine wants? see clip_dzdt() too.
end


"""
    eval_nn(density::NeuralDepthModel, z::FT, swe::FT, P::FT, T::FT, R::FT, qrel::FT, u::FT)::FT where {FT}
Helper function for evaluating the neural network in a pointwise manner over a `ClimaCore.Field`
and returning the output in a broadcastable way.
"""
function eval_nn(
    density::NeuralDepthModel,
    z::FT,
    swe::FT,
    P::FT,
    T::FT,
    R::FT,
    qrel::FT,
    u::FT,
    scf::FT,
)::FT where {FT}
    return density.z_model(
        SVector(z, swe_snow_area(swe, scf, z), qrel, R, u, T, P),
    )[1]
end

"""
    update_dzdt!(dzdt, density::NeuralDepthModel, Y, p, t)

Updates the dY.snow.Z field in places with the predicted change in snow depth (rate) given the model state `Y` and the `NeuralDepthModel`
density paramterization.
"""
function update_dzdt!(dzdt, density::NeuralDepthModel, Y, p)
    dzdt .=
        eval_nn.(
            Ref(density),
            Y.snow.Z,
            Y.snow.S,
            Y.snow.P_avg,
            Y.snow.T_avg,
            Y.snow.R_avg,
            Y.snow.Qrel_avg,
            Y.snow.u_avg,
            p.snow.snow_cover_fraction,
        )
end

"""
    clip_dZdt(S::FT, Z::FT, dSdt::FT, dZdt::FT, Δt::FT)::FT

A helper function which clips the tendency of Z such that
its behavior is consistent with that of S: if all snow melts
within a timestep, we clip the tendency of S so that it does
not become negative, and here we also clip the tendency of Z
so that depth does not become negative. Additionally, if the
tendencies of Z and S are such that we would encounter Z < S
(rho_snow > rho_liq, for Z and S over the snow-area, not
the grid cell area), we also clip the tendency.
"""
function clip_dZdt(
    S::FT,
    Z::FT,
    dSdt::FT,
    dZdt::FT,
    scf::FT,
    Δt::FT,
)::FT where {FT}
    min_scf = 0.05 #placeholder until Katherine and Andy discuss this
    #Case if S is set to zero:
    if (S + dSdt * Δt) <= eps(FT)
        return -Z / Δt
        #Case if Z would have been set to Z < S:
    elseif (Z + dZdt * Δt) * max(scf, min_scf) < (S + dSdt * Δt) #technically, the scf used here should be the scf of the following timestep, so how to handle this?
        #more stable form for very small Z, S
        return ((dSdt * Δt + S) / max(scf, min_scf) - Z) / Δt
    else
        return dZdt
    end
end

"""
    update_density_prog!(density::NeuralDepthModel, model::SnowModel, Y, p)

Updates all prognostic variables associated with density/depth given the current model state and the `NeuralDepthModel`
density paramterization.
"""
function update_density_prog!(
    density::NeuralDepthModel,
    model::SnowModel,
    dY,
    Y,
    p,
)
    update_dzdt!(dY.snow.Z, density, Y, p)

    # Now we clip the tendency so that Z stays within approximately physical bounds.
    @. dY.snow.Z = clip_dZdt(
        Y.snow.S,
        Y.snow.Z,
        dY.snow.S, #assumes dY.snow.S is updated (and clipped) before dY.snow.Z
        dY.snow.Z,
        p.snow.snow_cover_fraction,
        model.parameters.Δt,
    )

    @. dY.snow.P_avg = density.α * (abs(p.drivers.P_snow) - Y.snow.P_avg)
    @. dY.snow.T_avg =
        density.α *
        (p.drivers.T - model.parameters.earth_param_set.T_freeze - Y.snow.T_avg)
    @. dY.snow.R_avg = density.α * (p.drivers.SW_d - Y.snow.R_avg)
    dY.snow.Qrel_avg .=
        density.α .* (
            Thermodynamics.relative_humidity.(
                model.boundary_conditions.atmos.thermo_params,
                p.drivers.thermal_state,
            ) .- Y.snow.Qrel_avg
        )
    @. dY.snow.u_avg = density.α * (p.drivers.u - Y.snow.u_avg)
end
