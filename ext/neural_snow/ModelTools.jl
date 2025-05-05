module ModelTools
using Flux, LinearAlgebra
using DataFrames, Dates
export make_model,
    setoutscale!,
    setboundscale!,
    LinearModel,
    make_timeseries,
    trainmodel!,
    custom_loss,
    evaluate,
    paired_timeseries,
    save_predictive_model_weights,
    load_model_weights!,
    freeze_fixed_layers!,
    convert_model

include("./WebTools.jl")

""" 
    MulLayer {T<:AbstractMatrix}
Creates a custom Flux layer to reduce memory allocations and improve speed when no
bias nor activation term is needed, and only a matrix multiplication layer.
"""
struct MulLayer{T <: AbstractMatrix{<:AbstractFloat}}
    W::T
end


""" 
   (layer::MulLayer)(x::AbstractArray{<:AbstractFloat})
Defines the propogation behavior of the MulLayer and restricts the type
"""
function (layer::MulLayer)(x::AbstractArray{<:AbstractFloat})
    layer.W * x
end

#Sets up/endows the layer type with Flux layer properties (gradient tracking, etc.)
Flux.@layer MulLayer


"""
    connection_optimized(upperb, lowerb, pred, input)
Generator function to create the memory/speed-optimized connection-function for 
the neural network when provided with an upper-bound function `upperb(pred, input)`
and lower-bound function `lowerb(pred, input)`, where `pred` refers
to the set of predicted models by the predictive component of the model,
and `input` refers to the inputs of the model. The connection function is generated with a call
`connectfunc(pred, input) = connection_optimized(upperb, lowerb, pred, input)`
for usage within the model, see `make_model()` or `make_model_paper()`.
If one wishes to pass boundary functions taking additional arguments,
the below code must be modified to expand the function's arguments
and those passed to `upperb()` and `lowerb()`.
# Arguments
- `upperb::Function`: A function for calculating the upper-boundary value from `pred` and `input`
- `lowerb::Function`: A function for calculating the lower-boundary value from `pred` and `input`
- `pred::AbstractArray{<:AbstractFloat}`: The array of predictions from the predictive component.
- `input::AbstractArray{<:AbstractFloat}`: The array of inputs sent to the predictive component for each prediction.
    according with the predictive model structure, this array is N × K, where N is the number of features and K is the
    number of predictions.
"""
function connection_optimized(
    upperb::Function,
    lowerb::Function,
    pred::AbstractArray{T},
    input::AbstractArray{T},
)::AbstractArray{T} where {T <: AbstractFloat}
    return [
        upperb(pred, input)
        lowerb(pred, input)
        pred
    ]
end


"""
    paper_upper_bound(precip_idx, pred, input)
Generator function to create the memory/speed-optimized upper-bound function
from the paper when given which index of the input features pertains to the
precipitation value. This also serves as an example of one way to write a
boundary function for the network type introduced in the paper.
The upper-bound function is generated with a call
`upperb(pred, input) = paper_upper_bound(precip_idx, pred, input)`
within `make_model_paper()` given the index "precip_idx" (in the case of the paper
and our code defaults, this will be equal to 7).
The actual boundary value is 0 whenever there is no precipitation, otherwise,
it is the relu of the output prediction (0 if no precipitation, max(`p`, 0) otherwise
for prediction `p`).
Upper boundary functions can be specified without the need of a generator, for
example:
`upper_bound(pred, input) = input[1, :] .^2 .+ cos.(input[2, :]) ./ pred[3, :]` 
is a valid function for a neural network taking at least 2 inputs and outputting 
at least 3 values (where the third is nonzero), or
`upper_bound(pred, input) = myfunction1.(input[1, :]) .+ mean(mydata[:, "OFFSET"])`
is a valid function whenever `myfunction1()` and the `mydata` DataFrame are defined in
the same scope as the upper_bound (as well as any necessary libraries, like Statistics
for the `mean()` function)
`upper_bound(pred, input) = 0` is also a valid boundary as the Int is readily broadcastable
over the necessary output array size.
Remember to double-check before using macros like `@inbounds` or `@.` to further optimize functions
for speed. Row- vs Column-orientation matters for these functions (they should output a row).
If you would like to use a boundary function that makes use of additional 
dynamic/passed arguments (for example, additional variables from data/model outputs
that are not passed to the network inputs, or `flags` that alter the boundary
function, etc.), you must also modify the arguments and the call to `upper()`
within the `connection_optimized()` function to make sure these values are passed.
# Arguments
- `precip_idx::Int`: Specifies which row of the inputs provides the precipitation value (see `input`).
- `pred::AbstractArray{<:AbstractFloat}`: The array of predictions from the predictive component (1 × K for K samples).
- `input::AbstractArray{<:AbstractFloat}`: The array of inputs sent to the predictive component for each prediction.
    With the predictive model structure, this array is N × K, where N is the number of features and K is the
    number of predictions.
"""
paper_upper_bound(
    precip_idx::Int,
    pred::AbstractArray{<:AbstractFloat},
    input::AbstractArray{<:AbstractFloat},
) = @. (input[precip_idx, :]' > 0) * relu(pred)


"""
    paper_lower_bound(depth_idx, pred, input, i)
Generator function to create the memory/speed-optimized lower-bound function
from the paper when given which index of the input features pertains to the
depth (z) value. This also serves as an example of one way to write a
boundary function for the network type introduced in the paper.
The lower-bound function is generated with a call
`lowerb(pred, input) = paper_upper_bound(precip_idx, pred, input)`
within `make_model_paper()` given the index "z_idx" (in the case of the paper
and our code defaults, this will be equal to 1).
The returned boundary value is just the value of z, though the true lower bound
value is -z/Δt for model timestep Δt. For the paper model specifically, the 
scaling by -1 and by 1/Δt are handled elsewhere within the network for further speed/memory
improvements and to allow the user to set the timestep with the `settimescale!()`.
Lower boundary functions can be specified without the need of a generator, for
example:
`lower_bound(pred, input) = input[1, :] .^ 2 .+ cos.(input[2, :]) ./ pred[3, i]` 
is a valid function for a neural network taking at least 2 inputs and outputting 
at least 3 values (where the third is nonzero), or
`lower_bound(pred, input) = myfunction1.(input[1, :]) + mean.(mydata[:, "OFFSET"])`
is a valid function whenever `myfunction1()` and the `mydata` DataFrame are defined in
the same scope as the upper_bound (as well as any necessary libraries, like Statistics
for the `mean()` function)
`lower_bound(pred, input) = 0` is also a valid boundary as the Int is readily broadcastable
over the necessary output array size.
Remember to double-check before using macros like @inbounds or @. to further optimize functions
for speed. Row- vs Column-orientation matters for these functions (they should output a row).
If you would like to use a boundary function that makes use of additional 
dynamic/passed arguments (like additional variables from data/model outputs
that are not passed to the network inputs, or `flags` that alter the boundary
function, etc.), you must also modify the arguments and call to `lowerb()`
for the `connection_optimized()` function to make sure these values are passed.
# Arguments
- `depth_idx::Int`: Specifies which row of the inputs provides the depth value (see `input`).
- `pred::AbstractArray{<:AbstractFloat}`: The array of predictions from the predictive component.
- `input::AbstractArray{<:AbstractFloat}`: The array of inputs sent to the predictive component for each prediction.
    according with the predictive model structure, this array is N × K, where N is the number of features and K is the
    number of predictions.
"""
paper_lower_bound(
    depth_idx::Int,
    pred::AbstractArray{<:AbstractFloat},
    input::AbstractArray{<:AbstractFloat},
) = input[depth_idx, :]'


""" 
    convert_model(model, FT)
Provides a utility funciton to set a model's float type by passing
the type as an argument.
"""
function convert_model(model, FT::Type{<:AbstractFloat})
    if FT == Float32
        model = Flux.f32(model)
    elseif FT == Float64
        model = Flux.f64(model)
    else
        error(
            "Conversion of Flux Chain weights to the desired float-type is not supported. Please implement the desired conversion
      in ModelTools.convert_model() or a custom constructor for the desired type.",
        )
    end
    return model
end


"""
    make_model(pred_model, upperb, lowerb; ftype, in_scale, nfeatures)
Create a constrained predictive network that is bounded by the functions
`upperb` and `lowerb`. The function either takes `nfeatures` inputs or
specifies initial scaling weights for each feature `in_scales` so data can be
passed in raw/physically-viable form. The form of the internal predictive model
can be passed in the form of a Flux `Chain` argument.
The constraints are enforced with a fixed-weight layer that works for any values
output by `upperb` and `lowerb` as long as `upperb` ≥ `lowerb` for the same inputs
(e.g., `upperb` being (input)² and `lowerb` as √(input) are valid choices as long as the input
is never ≤ 1, or `upperb` = cos(input)+0.25 and `lowerb` = cos(input)-0.25 are also valid choices
over the entire domain even though `upperb` for some inputs is less than `lowerb` at other inputs).
If the boundary functions have additional properties like always being nonpositive or
nonnegative, the fixed weight layer can be reduced to further speed up computation - for
discussion on this see the paper, and for an exmample see the `NeuralDepthModel<:AbstractDensityModel` type in ClimaLand.Snow.

The user should specify either the number of features (`nfeatures`), or the scaling weights of each feature (`in_scale`) as an
argument for the model. If no scaling weights are specified, all model inputs will be scaled by a factor of 1. If both are provided,
the `in_scale` argument will take precedence (this argument allows for the fixing of input scalings so training and testing data can be input
in their original forms, and reduce the amount of preprocessing during model usage).

If building custom types with reduced fixed-weight layers or altnernative structures, the implementations of
functions like `setboundscale!()`, `setoutscale!()`, `freeze_fixed_layers!()`, or saving/loading utilities will likely require custom extension.

# Arguments
-  `pred_model`::Flux.Chain: The predictive chain composing the predictive model for this network
- `upperb::Function`: The upper bound function constraining the network output,
    see `paper_upper_bound() on how to construct`
- `lowerb::Function`: The lower bound function constraining the network output,
    see `paper_lower_bound() on how to construct`
- `ftype::Type`: Sets type of output model, currently `Float32` or `Float64` are supported.
- `nfeatures::Int`: indicates number of features.
- `in_scale::Vector{<:Real}`: Optional scaling constants for each input feature.
"""
function make_model(
    pred_model::Chain,
    upperb::Function,
    lowerb::Function;
    ftype::Type{FT} = Float32,
    in_scale::Union{Vector{<:Real}, Nothing} = nothing,
    nfeatures::Int = 0,
)::Chain where {FT <: AbstractFloat}
    if (nfeatures == 0) & isnothing(in_scale)
        error("Must specify either of `in_scale` or `nfeatures`\n")
    elseif (nfeatures != 0) & !isnothing(in_scale)
        @warn("Both `nfeatures` and `in_scale` are provided arguments to `make_model()`, defaulting to the provided value of `in_scale`\n")
        use_scales = Vector{ftype}(1 ./ in_scale)
    elseif !isnothing(in_scale)
        use_scales = Vector{ftype}(1 ./ in_scale)
    else
        use_scales = Vector{ftype}(ones(nfeatures))
    end
    #Explicitly set the index variables so any saved model does
    #not throw a reference error when loaded:
    let in_scales = use_scales,
        FT = ftype,
        upperbnd = upperb,
        lowerbnd = lowerb,
        predn = pred_model

        connectfunc(pred, input) =
            connection_optimized(upperbnd, lowerbnd, pred, input)

        # The general matrices for arbitary boundaries are used:
        # (see the paper on how to reduce these for certain boundary functions)
        # if making your own, the skip-connection outputs [upper, lower, p]ᵀ
        # and `setoutscale!()`, `setboundscale!()` might not work anymore, as
        # well as the saving/loading functions, all of which will likely require custom extensions.
        scale_and_relus = Matrix{FT}([1 0 0; -1 0 0; 0 1 0; 0 -1 0; 1 0 -1])
        get_min = Matrix{FT}([1 -1 -1 1 -1; 0 0 1 0 0; 0 0 0 1 0])
        final_mul = Matrix{FT}([1 1 -1])

        #The actual neural-network definition:
        model = Chain(
            pred = SkipConnection(
                Chain(; scale = x -> x .* in_scales, predn.layers...),
                connectfunc,
            ),
            apply_relus = Dense(scale_and_relus, false, relu),
            apply_upper = Dense(get_min, false, relu),
            sumlayer = MulLayer(final_mul),
        )
        return convert_model(model, FT)
    end
end


"""
    make_model_paper(nfeatures, n, depth_index, precipitation_index, ftype, in_scale)
Create the neural network structure from the paper, with initial scaling weights.
This does not set the predictive weight values - these must be loaded or the model can be
trained on additional data. This represents one specific class of the set of models that
could be made or represented with `make_model()`, though with additional optimization choices.

The model can be evaluated by passing a vector of input features, or a N × K matrix where
there are N features and K inputs to evaluate. This form also uses a reduced form of the
boundary enforcement layer for faster computation.

# Arguments
- `nfeatures::Int`: indicates number of features.
- `n::Int`: the value of the hyperparameter n.
- `depth_index::Int`: The index of the input feature vector pertaining to the depth (z) values,
    for our code the default value of this is 1.
- `precipitation_index::Int`: The index of the input feature vector pertaining to the
    precipitation values, for our code the default value of this is 7.
- `ftype::Type`: Sets type of utilized/output float. Default is `Float32`.
-  `in_scale::Vector{<:Real}`: Optional scaling constants for each input feature.
    Default values are set to those used in the paper.
"""
function make_model_paper(;
    nfeatures::Int = 7,
    n::Int = 4,
    depth_index::Int = 1,
    precipitation_index::Int = 7,
    ftype::Type{FT} = Float32,
    in_scale::Union{AbstractVector{FT}, Nothing} = FT.(
        [
            0.68659294 # z (m)
            0.25578135 # SWE (m)
            0.20772743 # relative humidity (0-1)
            76.2825 # solar radiation (W/m²)
            0.63011056 # wind speed (m/s)
            6.3486657 # air temp (⁰C)
            6.9992964f-8 # water-equiv rate of snowfall (m/s)
        ]
    ),
)::Chain where {FT <: AbstractFloat}
    if (precipitation_index > nfeatures) | (depth_index > nfeatures)
        error("The provided `nfeatures` is less than a 
        provided feature index, this will create out-of-bounds 
        errors during usage. Please revise either the depth or the 
        precipitation index, and make sure passed inputs are
        vectors of size `nfeatures`` or matrices of
        size (`nfeatures` × k) for k≥1.")
    end

    #Explicitly set the index variables so any saved model does
    #not throw a reference error when loaded:
    let nft = nfeatures,
        n_z = n,
        FT = ftype,
        p_idx = precipitation_index,
        z_idx = depth_index

        in_scales =
            (isnothing(in_scale)) ? Vector{FT}(ones(nfeatures)) :
            Vector{FT}(1 ./ in_scale)

        upperb(pred, input) = paper_upper_bound(p_idx, pred, input)
        lowerb(pred, input) = paper_lower_bound(z_idx, pred, input)
        connectfunc(pred, input) =
            connection_optimized(upperb, lowerb, pred, input)

        # These matrices set up the enforcing structure
        # for a nonpositive lower bound -z/Δt and a nonnegative
        # upper bound 1_{P>0}. The 1/Δt factor in the
        # lower bound will carried in the `scale_and_relus`
        # matrix and not the boundary function, and the -1 factor
        # is already built into the below matrix structure
        # to enable easier setting of the timescale and permit
        # a reduced structure. See the paper on how to set up
        # the enforcing structure for alternative bounds.
        scale_and_relus = Matrix{FT}([1 0 0; 0 1 0; 1 0 -1])
        get_min = Matrix{FT}([1 1 -1; 0 1 0])
        final_mul = Matrix{FT}([1 -1])

        #The actual neural-network definition:
        model = Chain(
            pred = SkipConnection(
                Chain(
                    scale = x -> x .* in_scales,
                    l1 = Dense(nft, n_z * nft, relu),
                    l2 = Dense(n_z * nft, nft, elu),
                    l3 = Dense(nft, 1),
                ),
                connectfunc,
            ),
            apply_relus = Dense(scale_and_relus, false, relu),
            apply_upper = Dense(get_min, false, relu),
            sumlayer = MulLayer(final_mul),
        )
        return convert_model(model, FT)
    end
end


"""
    setoutscale!(model, scale)

Set the output scaling parameter for model usage (i.e. rectifying scaling done on model input).
This may have unintended results or not work for models not created with `make_model()` or `make_model_paper()`.

Note: This should only be used for models made with `make_model()` or `make_model_paper()`.
This scaling constant is combined with a -1 factor that enforces the boundary for speed and memory purposes -
for this reason, we do not recommend first `extracting`/storing the weight this function refers to,
and then using that value with this function to `reset` the weight the intial/stored value, as this will have unintended
consequences without accounting for the extra factor. Only use this function to scale the output 
of the predictive component by the constant `scale`.

# Arguments
- `model::Chain`: the neural model to be scaled.
- `scale::Real`: the scaling parameter to return data to applicable units.
"""
function setoutscale!(model::Chain, scale::Real)
    FT = eltype(model[:apply_relus].weight)
    model[:apply_relus].weight[end, 3] = FT(-1 * scale)
end


"""
    setboundscale!(model, bound, scale)

Permits scaling of the upper or lower boundary function by the provided constant.
This may have unintended results or not work for models not created with `make_model()`.
For instance, in the neural snow depth paper (Charbonneau et. al., 2025), the lower bound
is -z/Δt, requiring a lower bound form of -z which is scaled by 1/Δt depending on the
utilized model timestep Δt.

Note: This function will only work for models made with `make_model()`, and will not work for those made with `make_model_paper()`, as reducing the form
of the matrices used to apply the bounds will cause uninteded results or errors with this function.
Users building custom models will likely need to write their own version of this function - e.g.,
the `settimescale!()` function for the `NeuralSnowDepth<:AbstractDensityModel` type or models made with `make_model_paper()`.
Like `setoutscale()!`, This scaling constant is combined with a -1 factor that enforces the boundary for speed and memory purposes -
for this reason, we do not recommend first `extracting`/storing the weights this function refers to,
and then using those value with this function to `reset` the weight the intial/stored value, as this will have unintended
consequences without accounting for the extra factor of -1. Only use this function to scale the bounds 
directly by the constant `scale`.

# Arguments
- `model::Chain`: the neural model to be scaled.
- `bound`::Symbol: specifies the upper or lower bound with :upper or :lower.
- `scale::Real`: the scaling parameter to return data to applicable units.
"""
function setboundscale!(model::Chain, bound::Symbol, scale::Real;)
    FT = eltype(model[:apply_relus].weight)
    if bound == :upper
        model[:apply_relus].weight[1, 1] = FT(1.0 * scale)
        model[:apply_relus].weight[2, 1] = FT(-1.0 * scale)
        model[:apply_relus].weight[5, 1] = FT(-1.0 * scale)
    elseif bound == :lower
        model[:apply_relus].weight[3, 2] = FT(1.0 * scale)
        model[:apply_relus].weight[4, 2] = FT(-1.0 * scale)
    else
        error("provided `bound` must be `:upper` or `:lower`.")
    end
end


"""
    settimescale!(model, dt)
Set the timescale parameter for model usage, for a model made with `make_model_paper()`.
Note: this should only be used on models made with `make_model_paper()`,
and not `make_model()`, or it will have unintended consequences. This function represents
a custom form of `setboundscale!()` to use with the custom model type made with `make_model_paper()`.

# Arguments
- `model::Chain`: the neural model to be used.
- `dt::Real`: the number of seconds per timestep for usage.
"""
function settimescale!(model::Chain, dt::Real)
    FT = eltype(model[:apply_relus].weight)
    model[:apply_relus].weight[2, 2] = FT(1.0 / dt)
end

"""
    freeze_fixed_layers!(model, opt_state)

Freezes the fixed layers of a model made with `make_model()` or `make_model_paper()` to prevent their alternation during training.
Required as of Flux version 0.15 instead of directly passing the Flux.params() of only the predictive component
of the model. The `model` argument is provided to enable easy extension of this method for other custom model
types, in order to make use of the `trainmodel!()` function for custom model types.

Note: this function will only work for models made with `make_model()` or `make_model_paper()`.

# Arguments
- `model::Chain`: the neural model utilized, made with `make_model()`.
- `opt_state`: The Flux optimization state for this model in a training pipeline.
"""
function freeze_fixed_layers!(model::Chain, opt_state)
    l = opt_state.layers
    @assert all(
        hasproperty.([model.layers], [:apply_relus, :apply_upper, :sumlayer]),
    ) "Provided model does not have the fixed layers for which freeze_fixed_layers!() was written, plase implement a custom extension for the supplied model"
    @assert all(hasproperty.([l], [:apply_relus, :apply_upper, :sumlayer])) "Provided Flux Setup State does not have the fixed layers for which freeze_fixed_layers!() was written, plase implement a custom extension for the supplied model"
    @assert hasproperty(model[:pred].layers.layers, :scale) "Provided Flux Setup State does not have the fixed input scaling vector for which this function was written, please implement a custom extension for the supplied model."
    @assert hasproperty(l[:pred].layers.layers, :scale) "Provided Flux Setup State does not have the fixed input scaling vector for which this function was written, please implement a custom extension for the supplied model."
    Flux.freeze!(l.apply_relus)
    Flux.freeze!(l.apply_upper)
    Flux.freeze!(l.sumlayer)
    Flux.freeze!(l.pred.layers.layers.scale)
end


"""
    LinearModel(data, vars, target; dtype, scale_const)

Create a linear regression model on a data frame for comparison to neural model.
Returns the coefficients for the model.

# Arguments
- `data::DataFrame`: The data set to be utilized.
- `vars::Vector{Symbol}`: The input variables to be used in the model creation.
-  `target::Symbol`: The target variable to be used in the model ceation.
- `dtype::Type`: Sets type, consistent with neural model. Default is `Float32`.
- `scale_const`: Optional scaling constant for model output. Default is 1.0.
"""
function LinearModel(
    data::DataFrame,
    vars::Vector{Symbol},
    target::Symbol;
    dtype::Type = Float32,
    scale_const = 1.0,
)
    X = Matrix{dtype}(select(data, vars))
    y = Vector{dtype}(data[!, target]) ./ dtype(scale_const)
    constants = [X ones(nrow(data))] \ y
    return Vector{dtype}(constants .* scale_const)
end


"""
    LinearModel(x_train, y_train; dtype, scale_const)

Create a linear regression model on a training matrix for comparison to neural model.
Returns the coefficients for the model.
**Note: using the same matrices input to the neural model will require a transpose of `x_train`, `y_train`

# Arguments
- `x_train::Matrix`: The input to be utilized.
- `y_train::Vector`: The output data to be utilized.
- `dtype::Type`: Sets type, consistent with neural model. Default is `Float32`.
- `scale_const`: Optional scaling constant for model output. Default is 1.0.
"""
function LinearModel(
    x_train::Matrix,
    y_train::Vector;
    dtype::Type = Float32,
    scale_const = 1.0,
)
    #using x_train from neural input will require a transpose of x_train, y_train
    return Vector{dtype}(
        ([x_train ones(size(x_train)[1])] \ y_train) .* scale_const,
    )
end


"""
    evaluate(model, input)

Evaluate a created model on a given input vector.

# Arguments
- `model::Chain`: A neural model to be used for prediction.
- `input`: The input data used to generate a prediction.
"""
function evaluate(model::Chain, input)
    return model(input)
end


"""
    evaluate(model, input)

Evaluate a created model on a given input vector.

# Arguments
- `model::Vector{<:Real}`: Linear regression coefficients used for prediction.
- `input`: The input data used to generate a prediction.
"""
function evaluate(model::Vector{<:Real}, input)
    #requires input matrix to be the same orientation as that for the neural model
    return model[1:(end - 1)]' * input .+ model[end]
end


"""
    make_timeseries(model, timeseries, dt; predictvar, timevar, inputvars, dtype, hole_thresh)

Generate a predicted timeseries given forcing data and the timestep present in that data (holes acceptable).
This function can be used for any predictive model regardless of its structure as long as the function has a defined
`evaluate(model, input)` method.

# Arguments
- `model`: The model used for forecasting (can be any model with a defined `evaluate` call).
- `timeseries::DataFrame`: The input data frame used to generate predictions, including a time variable.
- `dt::Period`: The unit timestep present in the dataframe (i.e. daily dataframe, `dt = Day(1)` or `Second(86400)`).
- `predictvar::Symbol`: The variable to predict from the timeseries. Default is `:z`.
- `timevar::Symbol`: The variable giving the time of each forcing. Default is `:date`.
- `inputvars::Vector{Symbol}`: The variables (in order), to extract from the data to use for predictions.
Default is `[:z, :SWE, :rel_hum_avg, :sol_rad_avg, :wind_speed_avg, :air_temp_avg, :dprecipdt_snow]` like the paper.
- `dtype::Type`: The data type required for input to the model. Default is `Float32`.
- `hole_thresh::Int`: The acceptable number of "holes" in the timeseries for the model to skip over. Default is 30.
"""
function make_timeseries(
    model,
    timeseries::DataFrame,
    dt::Period;
    predictvar::Symbol = :z,
    timevar::Symbol = :date,
    inputvars::Vector{Symbol} = [
        :z,
        :SWE,
        :rel_hum_avg,
        :sol_rad_avg,
        :wind_speed_avg,
        :air_temp_avg,
        :dprecipdt_snow,
    ],
    dtype::Type = Float32,
    hole_thresh::Int = 30,
)
    forcings = Matrix{dtype}(select(timeseries, inputvars))
    pred_idx = findfirst(inputvars .== predictvar)
    check_dates =
        (timeseries[2:end, timevar] - timeseries[1:(end - 1), timevar]) ./ dt
    pred_vals = zeros(nrow(timeseries) - 1)
    pred_series = zeros(nrow(timeseries))
    pred_series[1] = timeseries[1, predictvar]
    countresets = 0
    for j in 2:length(pred_series)
        input = forcings[j - 1, :]
        input[pred_idx] = pred_series[j - 1]
        pred = evaluate(model, input)[1]
        pred_vals[j - 1] = pred
        nperiods = check_dates[j - 1]
        new_val = pred_series[j - 1] + nperiods * Dates.value(Second(dt)) * pred
        pred_series[j] =
            (nperiods <= hole_thresh) ? max(0.0, new_val) :
            timeseries[j, predictvar]  #the "max" is only in the case of holes for a non-negative system, and does not generalize
        #pred_series[j] = (nperiods <= hole_thresh) ? new_val : true_series[j]  #for showing no thresholds
        if nperiods > hole_thresh
            countresets += 1
        end
    end
    return pred_series, pred_vals, countresets
end


"""
    paired_timeseries(zmodel, swemodel, timeseries, dt; timevar, z_inputvars, SWE_inputvars, dtype, hole_thresh)

Generate a predicted timeseries given forcing data and the timestep present in that data (holes acceptable), with a coupled
z and SWE model. Used in the paper for combined timeseries.

# Arguments
- `zmodel`: The model used for depth forecasting (can be any model with a defined `evaluate` call).
- `swemodel`: The model used for SWE forecasting (can be any model with a defined `evaluate` call).
- `timeseries::DataFrame`: The input data frame used to generate predictions, including a time variable.
- `dt::Period`: The unit timestep present in the dataframe (i.e. daily dataframe, `dt = Day(1)` or `Second(86400)`).
- `timevar::Symbol`: The variable giving the time of each forcing. Default is `:date`.
- `z_inputvars::Vector{Symbol}`: The variables (in order), to extract from the data to use for depth predictions.
Default is `[:z, :SWE, :rel_hum_avg, :sol_rad_avg, :wind_speed_avg, :air_temp_avg, :dprecipdt_snow]` like the paper.
- `SWE_inputvars::Vector{Symbol}`: The variables (in order), to extract from the data to use for SWE predictions.
Default is `[:z, :SWE, :rel_hum_avg, :sol_rad_avg, :wind_speed_avg, :air_temp_avg, :dprecipdt_snow]` like the paper.
- `dtype::Type`: The data type required for input to the model. Default is `Float32`.
- `hole_thresh::Int`: The acceptable number of "holes" in the timeseries for the model to skip over. Default is 30.
"""
function paired_timeseries(
    zmodel,
    swemodel,
    timeseries::DataFrame,
    dt::Period;
    timevar::Symbol = :date,
    z_inputvars::Vector{Symbol} = [
        :z,
        :SWE,
        :rel_hum_avg,
        :sol_rad_avg,
        :wind_speed_avg,
        :air_temp_avg,
        :dprecipdt_snow,
    ],
    SWE_inputvars::Vector{Symbol} = [
        :z,
        :SWE,
        :rel_hum_avg,
        :sol_rad_avg,
        :wind_speed_avg,
        :air_temp_avg,
        :dprecipdt_snow,
    ],
    dtype::Type = Float32,
    hole_thresh::Int = 30,
)
    zforcings = Matrix{dtype}(select(timeseries, z_inputvars))
    sweforcings = Matrix{dtype}(select(timeseries, SWE_inputvars))
    z_zidx = findfirst(z_inputvars .== :z)
    z_sweidx = findfirst(z_inputvars .== :SWE)
    swe_zidx = findfirst(SWE_inputvars .== :z)
    swe_sweidx = findfirst(SWE_inputvars .== :SWE)
    check_dates =
        (timeseries[2:end, timevar] - timeseries[1:(end - 1), timevar]) ./ dt
    pred_dzdts = zeros(nrow(timeseries) - 1)
    pred_dswedts = zeros(nrow(timeseries) - 1)
    pred_zs = zeros(nrow(timeseries))
    pred_swes = zeros(nrow(timeseries))
    pred_zs[1] = timeseries[1, :z]
    pred_swes[1] = timeseries[1, :SWE]
    countresets = 0
    for j in 2:length(pred_zs)
        zinput = zforcings[j - 1, :]
        sweinput = sweforcings[j - 1, :]
        zinput[z_zidx] = pred_zs[j - 1]
        zinput[z_sweidx] = pred_swes[j - 1]
        sweinput[swe_zidx] = pred_zs[j - 1]
        sweinput[swe_sweidx] = pred_swes[j - 1]
        dzdtpred = evaluate(zmodel, zinput)[1]
        dswedtpred = evaluate(swemodel, sweinput)[1]
        pred_dzdts[j - 1] = dzdtpred
        pred_dswedts[j - 1] = dswedtpred
        nperiods = check_dates[j - 1]
        new_z = pred_zs[j - 1] + nperiods * Dates.value(Second(dt)) * dzdtpred
        new_swe =
            pred_swes[j - 1] + nperiods * Dates.value(Second(dt)) * dswedtpred
        pred_swes[j] =
            (nperiods <= hole_thresh) ? max(0.0, new_swe) : timeseries[j, :SWE]
        pred_zs[j] =
            (nperiods <= hole_thresh) ? max(pred_swes[j], new_z) :
            timeseries[j, :z]
        if nperiods > hole_thresh
            countresets += 1
        end
    end
    return pred_zs, pred_swes, pred_dzdts, pred_dswedts, countresets
end


"""
    custom_loss(x, y, model, n1, n2)

Creates a loss function to be used during training specified by two hyperparameters n1, n2, as outlined in the paper.

# Arguments
- `x`: the input values.
- `y`: output values to compare to.
- `model`: the model over which to evaluate the loss function.
- `n1::Int`: the hyperparameter dictating the scaling of mismatch error.
- `n2::Int`: the hyperparameter dictating the weighting of a given mismatch error by the target magnitude.
"""
function custom_loss(x, y, model, n1, n2)
    sum(abs.((model(x) .- y)) .^ n1 .* (1 .+ abs.(y) .^ n2)) ./ length(y)
end


"""
    trainmodel!(model, x_train, y_train, n1, n2; nepochs, opt, verbose, cb)

A training function for a neural model, permitting usage of a callback function.
If a custom model is being used, make sure the model type has an appropriately defined/extended `freeze_fixed_layers!()` function,
as of Flux version 0.15.

# Arguments
- `model`: the model used for training.
- `x_train`: the input training data to be used in training.
- `y_train`: the target data to be used in training.
- `n1`: the scaling hypermarameter used to generate custom loss functions.
- `n2`: the weighting hyperparameter used to generate custom loss functions.
- `nepochs::Int`: the number of epochs. Default is 100.
- `nbatch::Int`: The number of data points to be used per batch. Default is 64.
- `opt`: the Flux optimizer to be used. Default is `RMSProp()`.
- `verbose::Bool`: indicates whether to print the training loss every 10 epochs
- `cb`: Allows utlization of a callback function (must take no required
input arguments, but default optional args are permitted). Default is Nothing.
"""
function trainmodel!(
    model,
    x_train,
    y_train,
    n1,
    n2;
    nepochs::Int = 100,
    nbatch = 64,
    opt = Flux.Optimisers.RMSProp(),
    verbose = false,
    cb = Nothing,
)
    train_loader = Flux.DataLoader(
        (x_train, y_train),
        batchsize = nbatch,
        partial = false,
        shuffle = true,
    )
    loss(_model, x, y) = custom_loss(x, y, _model, n1, n2)
    opt_state = Flux.setup(opt, model)
    freeze_fixed_layers!(model, opt_state)
    for epoch in 1:nepochs
        for (x, y) in train_loader
            Flux.train!(loss, model, [(x, y)], opt_state)
        end
        if verbose & (epoch % 10 == 0)
            print(
                "Epoch: ",
                epoch,
                " | training loss: ",
                loss(model, x_train, y_train),
                "\n",
            )
        end
        cb()
    end
end


"""
   save_predictive_model_weights(path, model)

A helper/utility function for writing the weights of the PREDICTIVE COMPONENT
of a model made with `make_model()` or `make_model_paper()` to a human-parseable text file.

This only saves the weights of the predictive component of the network, and only
if the predictive component is a chain of named `Dense` layers, in the same manner
as the paper by Charbonneau et. al. (2025). For custom models or alternative constructions
we refer the reader to the Flux standard Saving and Loading page: https://fluxml.ai/Flux.jl/stable/guide/saving/
This provides parameter files that can be readily used with other coding languages and is human readable,
as well as working independent of the Flux package version.

# Arguments
- `path`::String: The filepath for which to save the text file.
- `model`: the model for which to save the predictive parameters.
"""
function save_predictive_model_weights(path::String, model::Chain)
    layers = model[:pred].layers.layers
    scalings = layers.scale.in_scales
    open(path, "w") do file
        write(file, "SCALING:\n")
        println(file, join(scalings, " "))
        for layer in collect(keys(layers))[2:end]
            write(file, "\n" * string(layer) * "-WEIGHT:\n")
            for row in eachrow(layers[layer].weight)
                println(file, join(row, " "))
            end
            write(file, "\n" * string(layer) * "-BIAS:\n")
            println(file, join(layers[layer].bias, " "))
        end
        write(file, "\nFINALSCALE:\n")
        println(file, -model[:apply_relus].weight[end, 3])
    end
end


"""
   load_model_weights!(path, model)

A helper/utility function for loading the weights of the PREDICTIVE COMPONENT
into a model made with `make_model()` or `make_model_paper()` from a text file (or web URL) for a model
with the same structure made with `save_predictive_model_weights()`.

Note this will only work for models matching the structure as those defined
in the text file, and is written specifically for models made with models
generated with `make_model()` or `make_model_paper()`. Calling this function will not work for custom
model types, or for models having different structures than those which had their
parameters saved with `save_predictive_model_weights()`. For custom models or alternative constructions
we refer the reader to the Flux standard Saving and Loading page: https://fluxml.ai/Flux.jl/stable/guide/saving/
This provides a loading utility that can be readily used with weight files generated from other coding languages,
and in a way that is stable to changing versions of the Flux package.

# Arguments
- `path`::String: The filepath (or URL) from where to get the text file.
- `model`: the model for which to load in the predictive parameters.
- `FT` : Optional arg for passing the right Float type if not inferred from passed model.
"""
function load_model_weights!(filepath::String, model::Chain; FT = Float32)
    T =
        hasproperty(model.layers, :apply_relus) ?
        eltype(model[:apply_relus].weight) : FT
    data = Dict()  # Dictionary to store the parsed data
    # Check if the filepath is a URL
    if startswith(filepath, "http://") || startswith(filepath, "https://")
        file_io = IOBuffer(read_webpage_body(filepath))
    else
        file_io = open(filepath, "r")  # Open as a regular file
    end
    try
        current_name = nothing
        current_values = []  # Temporary storage for current block

        for line in readlines(file_io)
            line = strip(line)
            if isempty(line) && current_name !== nothing
                # If an empty line is encountered, finalize the current item
                if length(current_values[1]) == 1 #one column of data (vector)
                    data[current_name] = vcat(current_values...)
                else # It's a matrix (multiple rows)
                    data[current_name] = hcat(current_values...)'
                end
                current_name = nothing
                current_values = []
            elseif isempty(line)
                continue
            elseif isnothing(current_name)
                current_name = line
            else
                push!(current_values, parse.(T, split(line)))
            end
        end
        # Handle the last item (if file doesn't end with a blank line)
        if current_name !== nothing
            if length(current_values[1]) == 1
                data[current_name] = vcat(current_values...)
            else
                data[current_name] = hcat(current_values...)'
            end
        end
    finally
        # Close the file only if it's an actual file, not an IOBuffer
        if isa(file_io, IO)
            close(file_io)
        end
    end
    layers = model[:pred].layers.layers
    layers[:scale].in_scales .= data["SCALING:"]'
    for layer in collect(keys(layers))[2:end]
        layers[layer].weight .= data[string(layer) * "-WEIGHT:"]
        layers[layer].bias .= data[string(layer) * "-BIAS:"]'
    end
    model[:apply_relus].weight[end, 3] = -1 * data["FINALSCALE:"][1]
end

end
