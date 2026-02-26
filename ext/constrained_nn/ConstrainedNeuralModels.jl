module ConstrainedNeuralModels

using Flux, JLD2, Adapt #weak dependencies for the extension (plus InteractiveUtils for `methodswith()`)
using StaticArrays #ClimaLand dependency
using ClimaLand #only here to be able to grab its dependencies and version number

export ConstrainedNeuralModel,
    set_predictive_model_out_scale!,
    scale_model!,
    save_model,
    load_model,
    make_static_model,
    make_dynamic_model,
    trainmodel!,
    inspect_model_metadata,
    build_model_bound_documentation,
    build_model_API,
    convert_model,
    @bound,
    @bound_type

"""
    InputWeighting{FT <: AbstractFloat}

Defines a type for modes of input weighting to a ConstrainedNeuralModel - that is,
whether inputs will require scaling or not.
"""
abstract type InputWeighting{FT <: AbstractFloat} end

"""
    NoScaling{FT} <: InputWeighting{FT}

Establishes the style of input weighting where no scaling is needed.
"""
struct NoScaling{FT <: AbstractFloat} <: InputWeighting{FT} end

"""
    ConstScaling{FT <: AbstractFloat, T <: AbstractArray{FT}} <: InputWeighting{FT}

Establishes the style of input weighting where each input is scaled by a specified constant.
"""
struct ConstScaling{FT <: AbstractFloat, T <: AbstractArray{FT}} <:
       InputWeighting{FT}
    in_scales::T
end

"""
    ConstraintType

Defines a type for the different constraint possibilities of a ConstrainedNeuralModel
(an upper bound, a lower bound, or both)
"""
abstract type ConstraintType end

"""
    UpperOnly{F} <: ConstraintType

Establishes the style of one-sided constraint serving to threshold the upper
limit of a model output.
"""
struct UpperOnly{F} <: ConstraintType
    bound::F
end

"""
    output_dim(::UpperOnly)::Int = 2

Sets the dimension output size for an UpperOnly constraint (for model output construction).
"""
@inline output_dim(::UpperOnly)::Int = 2

#= Leave out as we won't automate checks anymore; using the default constructor
# adding this back in would require "UpperOnly" -> "buildUpper" in `buildConstraints()`
"""
    buildUpper(bound)

Constructs an UpperOnly constraint from a provided bound.
Passed bounds should be defined with the `@bound` macro on their bound definition,
to ensure the bound is compliant for usage within ClimaLand.
"""
function buildUpper(bound)::UpperOnly
    if is_valid_bound(bound)
        return UpperOnly{typeof(bound)}(bound)
    else
        error(
            "Valid upper bound functions must be specified with @bound macro.",
        )
    end
end
=#

"""
    LowerOnly{F} <: ConstraintType

Establishes the style of one-sided constraint serving to threshold the lower
limit of a model output.
"""
struct LowerOnly{F} <: ConstraintType
    bound::F
end

"""
    output_dim(::LowerOnly)::Int = 2

Sets the dimension output size for an LowerOnly constraint (for model output construction).
"""
@inline output_dim(::LowerOnly)::Int = 2

#= Leave out as we won't automate checks anymore; using the default constructor
# adding this back in would require "LowerOnly" -> "buildLower" in `buildConstraints()`
"""
    buildLower(bound)

Constructs an LowerOnly constraint from a provided bound.
Passed bounds should be defined with the `@bound` macro on their bound definition,
to ensure the bound is compliant for usage within ClimaLand.
"""
function buildLower(bound)::LowerOnly
    if is_valid_bound(bound)
        return LowerOnly{typeof(bound)}(bound)
    else
        error("Valid lower bound function must be specified with @bound macro.")
    end
end
=#

"""
    TwoSided{F1, F2} <: ConstraintType

Establishes the style of a two-sided constraint, serving to threshold the upper
and lower limits of a model output.
"""
struct TwoSided{F1, F2} <: ConstraintType
    upper_bound::F1
    lower_bound::F2
end

"""
    output_dim(::TwoSided)::Int = 3

Sets the dimension output size for an TwoSided constraint (for model output construction).
"""
@inline output_dim(::TwoSided)::Int = 3

#= Leave out as we won't automate checks anymore; using the default constructor
# adding this back in would require "TwoSided" -> "buildTwoSided" in `buildConstraints()`
"""
    buildTwoSided(bound)

Constructs a TwoSided constraint from provided bounds.
Passed bounds should be defined with the `@bound` macro on their bound definitions,
to ensure the bounds are compliant for usage within ClimaLand.
"""
function buildTwoSided(upper_bound, lower_bound)::TwoSided
    if is_valid_bound(upper_bound) && is_valid_bound(lower_bound)
        return TwoSided{typeof(upper_bound), typeof(lower_bound)}(
            upper_bound,
            lower_bound,
        )
    else
        error("Valid bound functions must be specified with @bound macro")
    end
end
=#

"""
    ScaleOutput{FT<:AbstractFloat, T <: AbstractVector{FT}}

Establishes a custom type to be used for the output scaling of a ConstrainedNeuralModel.
Consists of a 1-element Vector, in order to preserve allocation-less multiplication
when working with a fixed model.
"""
struct ScaleOutput{FT <: AbstractFloat, T <: AbstractVector{FT}}
    sc::T
end

"""
    ScaleOutput(x::FT) where {FT <: AbstractFloat}

Constructs a ScaleOutput layer from the provided scaling constant.
"""
function ScaleOutput(x::FT)::ScaleOutput{FT} where {FT <: AbstractFloat}
    sc = [x]
    return ScaleOutput{FT, typeof(sc)}(sc)
end

"""
    function _get_scale(x::ScaleOutput{FT})::FT where {FT<:AbstractFloat}

Internal function for returning the ScaleOuput's scaling constant.
"""
@inline function _get_scale(x::ScaleOutput{FT})::FT where {FT <: AbstractFloat}
    return x.sc[1]
end

"""
    (layer::ScaleOutput{FT})(x::AbstractArray{FT})::AbstractArray{FT} where {FT <: AbstractFloat}

Establishes the behavior of a ScaleOuput layer for generic inputs.
"""
@inline function (layer::ScaleOutput{FT})(
    x::T,
)::T where {FT <: AbstractFloat, T <: AbstractArray{FT}}
    return _get_scale(layer) .* x
end

#Establishes the ScaleOuput type as a valid Flux layer:
Flux.@layer ScaleOutput

"""
    MulLayer{FT<:AbstractFloat, T <: AbstractMatrix{FT}}

Defines a Flux model layer (a functor) for weight multiplication when no bias term is
necessary, reducing the computational resources compared to a standard Flux Dense layer.
"""
struct MulLayer{FT <: AbstractFloat, T <: AbstractMatrix{FT}}
    W::T
end

"""
    (layer::MulLayer{FT})(x::AbstractArray{FT})::AbstractArray{FT} where {FT<:AbstractFloat}

Defines MulLayer behavior for a generic input.
"""
@inline function (layer::MulLayer{FT})(
    x::T,
) where {FT <: AbstractFloat, T <: AbstractArray{FT}}
    return layer.W * x
end

#Establishes the MulLayer as a Flux layer now that its behavior is defined:
Flux.@layer MulLayer

"""
    get_bounds(c::ConstraintType, pred, input)

Returns the output of a ConstraintType for a given model prediction and input, an internal utility
in defining the information flow of an input through a ConstrainedNeuralModel.
Note: if pred, input are not passed, this function is additionally defined to just return the bounds
themselves, not their evaluations for a given prediction and input.
"""
@inline get_bounds(c::Union{UpperOnly, LowerOnly}, pred, input) =
    (c.bound(pred, input),)
@inline get_bounds(c::TwoSided, pred, input) =
    (c.upper_bound(pred, input), c.lower_bound(pred, input))

"""
    boundary_connection(
        constraints::ConstraintType,
        pred::VecOrMat{FT},
        input::AbstractArray{FT},
    )::VecOrMat{FT} where {FT <: AbstractFloat}

Defines the behavior connecting the output of a predictive model to the fixed layers which apply the constraints.
Analagous to defining a custom output function for a Flux SkipConnection, instead of the default "vcat" or "+" operations.
This generic version works for many possible input types for the user, defined for prototyping on CPU.
"""
@inline function boundary_connection(
    constraints::ConstraintType,
    pred::Union{FT, <:VecOrMat{FT}},
    input::AbstractArray{FT},
)::VecOrMat{FT} where {FT <: AbstractFloat}
    return [
        get_bounds(constraints, pred, input)...
        pred
    ]
end

"""
    boundary_connection(
        constraints::ConstraintType,
        pred::SMatrix{1, N, FT, N},
        input::StaticArray{S, FT},
    )::SMatrix{output_dim(constraints), N, FT} where {S, N, FT <: AbstractFloat}

Defines the behavior connecting the output of a predictive model to the fixed layers which apply the constraints.
Analagous to defining a custom output function for a Flux SkipConnection, instead of the default "vcat" or "+" operations.
This specialized version works optimizing performance for batched inputs on GPU. It specifically requires boundary functions
to output a 1-dimensional SMatrix.
"""
@inline function boundary_connection(
    constraints::ConstraintType,
    pred::SMatrix{1, N, FT, N},
    input::StaticArray{S, FT},
)::SMatrix{output_dim(constraints), N, FT} where {S, N, FT <: AbstractFloat}
    SMatrix{output_dim(constraints), N, FT}(
        vcat(get_bounds(constraints, pred, input)..., pred),
    )
end

"""
    boundary_connection(
        constraints::ConstraintType,
        pred::SVector{1, FT},
        input::SArray{S, FT},
    )::SMatrix{output_dim(constraints), 1, FT} where {S, FT <: AbstractFloat}

Defines the behavior connecting the output of a predictive model to the fixed layers which apply the constraints.
Analagous to defining a custom output function for a Flux SkipConnection, instead of the default "vcat" or "+" operations.
This specialized version works optimizing performance for single inputs on GPU. It specifically requires boundary functions
to output a single float.
"""
@inline function boundary_connection(
    constraints::ConstraintType,
    pred::SVector{1, FT},
    input::SArray{S, FT},
)::SMatrix{output_dim(constraints), 1, FT} where {S, FT <: AbstractFloat}
    return SMatrix{output_dim(constraints), 1, FT}(
        get_bounds(constraints, pred, input)...,
        pred[1],
    )
end

"""
    convert_model(model, FT::Type{<:AbstractFloat})

Utility function for converting models or structures to be of a consistent float type. Depending
on your custom types, this function might not work on an entire ConstrainedNeuralModel.
"""
function convert_model(model, FT::Type{<:AbstractFloat})
    if FT == Float32
        return fmap(f32, model)
    elseif FT == Float64
        return fmap(f64, model)
    else
        throw(
            ArgumentError(
                """
    convert_model(): Conversion of parameters to the desired float-type is
    not supported. Please implement the desired conversion in ConstrainedNeuralModels.convert_model()
    or a custom constructor for the desired type.
    """,
            ),
        )
    end
end

"""
    ConstrainedNeuralModel{
        FT <: AbstractFloat,
        S <: InputWeighting{FT},
        TP <: Val,
        C <: ConstraintType,
        NN1,
        NN2 <: Chain,
        IL <: AbstractMatrix{FT},
        OS <: ScaleOutput{FT}
    }

Defines the ConstrainedNeuralModel type, taking an InputWeighting and a ConstraintType, as well as a predictive model
as well as a Flux Chain to act as the fixed layers that apply the constraints. An output scaling
is used to aid scaling for model training, as well as an indication if the user is using custom fixed layers compared
to the default, as well as whether the constraint parameterizations are trainable. The initial fixed layer weight
is likewise stored to enable resetting after rescaling.
"""
struct ConstrainedNeuralModel{
    FT <: AbstractFloat,
    S <: InputWeighting{FT},
    TP <: Val,
    C <: ConstraintType,
    NN1,
    NN2 <: Chain,
    IL <: AbstractMatrix{FT},
    OS <: ScaleOutput{FT},
}
    predictive_model::NN1
    constraints::C
    scaling::S
    out_scale::OS
    fixed_layers::NN2
    initial_fixed_layer::IL
    using_default_fixed_layers::Bool
    trainable_constraints::TP
end

"""
    (model::ConstrainedNeuralModel{<:ConstScaling})(x::AbstractArray{<:AbstractFloat})::AbstractArray{<:AbstractFloat}

Defines the behavior of a ConstrainedNeuralModel type for an input, when input scaling is required. Inputs are scaled,
and the output of the boundary functions and predictive model are combined and passed through the fixed layer which
applies the constraints.
"""
function (model::ConstrainedNeuralModel{FT, C})(
    x,
) where {FT <: AbstractFloat, C <: ConstScaling{FT}}
    return model.fixed_layers(
        boundary_connection(
            model.constraints,
            model.out_scale(
                model.predictive_model(x .* model.scaling.in_scales),
            ),
            x,
        ),
    )
end

"""
    (model::ConstrainedNeuralModel{<:NoScaling})(x::AbstractArray{<:AbstractFloat})::AbstractArray{<:AbstractFloat}

Defines the behavior of a ConstrainedNeuralModel type for an input, when no scaling is required. This serves to
reduce computational cost compared to rescaling inputs by a factor of 1. The output of the boundary functions
and predictive model are combined and passed through the fixed layer which applies the constraints.
"""
function (model::ConstrainedNeuralModel{FT, C})(
    x,
) where {FT <: AbstractFloat, C <: NoScaling{FT}}
    return model.fixed_layers(
        boundary_connection(
            model.constraints,
            model.out_scale(model.predictive_model(x)),
            x,
        ),
    )
end

#Establishes the ConstrainedNeuralModel as a Flux layer now that its behavior is defined:
Flux.@layer ConstrainedNeuralModel

"""
    Flux.trainable(m::ConstrainedNeuralModel{S, Val{false}}) where {S<:InputWeighting}

Sets the trainable parameters of the ConstrainedNeuralModel to just the predictive model,
when constraints are indicated as untrainable, to avoid the need to freeze layers prior to training.
"""
Flux.trainable(
    m::ConstrainedNeuralModel{FT, S, Val{false}},
) where {FT <: AbstractFloat, S <: InputWeighting{FT}} = (; m.predictive_model)

"""
    Flux.trainable(m::ConstrainedNeuralModel{S, Val{true}, C}) where {S<:InputWeighting, C<:Union{<:UpperOnly, <:LowerOnly}}

Sets the trainable parameters of the ConstrainedNeuralModel to just the predictive model and the single bound,
when constraints are indicated as trainable, to avoid the need to freeze/thaw layers prior to training.
"""
Flux.trainable(
    m::ConstrainedNeuralModel{FT, S, Val{true}, C},
) where {
    FT <: AbstractFloat,
    S <: InputWeighting{FT},
    C <: Union{<:UpperOnly, <:LowerOnly},
} = (; m.predictive_model, m.constraints.bound)

"""
    Flux.trainable(m::ConstrainedNeuralModel{S, Val{true}, C}) where {S<:InputWeighting, C<:TwoSided}

Sets the trainable parameters of the ConstrainedNeuralModel to just the predictive model and the single bound,
when constraints are indicated as trainable, to avoid the need to freeze/thaw layers prior to training.
"""
Flux.trainable(
    m::ConstrainedNeuralModel{FT, S, Val{true}, C},
) where {FT <: AbstractFloat, S <: InputWeighting{FT}, C <: TwoSided} =
    (; m.predictive_model, m.constraints.upper_bound, m.constraints.lower_bound)

"""
    buildConstraints(upper_bound, lower_bound)::ConstraintType

Creates a ConstraintType from passed bounds, ensuring bounds exist to create a ConstrainedNeuralModel.
"""
#If testing bounds is made required again -> change constructors to equivalent "build__()" functions.
function buildConstraints(upper_bound, lower_bound)::ConstraintType
    if isnothing(upper_bound) && isnothing(lower_bound)
        throw(
            ArgumentError(
                "You must supply at least one constraining function to a ConstrainedNeuralModel.",
            ),
        )
    elseif isnothing(upper_bound) && !isnothing(lower_bound)
        return LowerOnly(lower_bound)
    elseif !isnothing(upper_bound) && isnothing(lower_bound)
        return UpperOnly(upper_bound)
    else
        return TwoSided(upper_bound, lower_bound)
    end
end

"""
    default_fixed_layers(constraint::LowerOnly, FT::Type{<:AbstractFloat})::Chain

Creates the default fixed-weight layers for applying a lower bound to an output. The general matrices for
arbitrary boundaries are used on output [f, p]ᵀ, where f is the bound function output and p is the output
of the predictive model. See the NeuralDepthModel paper (https://doi.org/10.1175/AIES-D-24-0040.1) for
discussion on how to reduce these for particular boundary functions for reduced computational cost.
"""
function default_fixed_layers(
    constraint::LowerOnly,
    FT::Type{<:AbstractFloat},
)::Chain
    get_relus = Matrix{FT}([1 0; -1 0; -1 1]) # (relu(f), relu(-f), relu(p - f))
    final_mul = Matrix{FT}([1 -1 1])          # outputs relu(f) - relu(-f) + relu(p-f) = f + relu(p - f) = max(f, p)
    return Chain(
        apply_relus = Dense(get_relus, false, relu),
        sumlayer = MulLayer(final_mul),
    )
end

"""
    default_fixed_layers(constraint::UpperOnly, FT::Type{<:AbstractFloat})::Chain

Creates the default fixed-weight layers for applying an upper bound to an output. The general matrices for
arbitrary boundaries are used on output [f, p]ᵀ, where f is the bound function output and p is the output
of the predictive model. See the NeuralDepthModel paper (https://doi.org/10.1175/AIES-D-24-0040.1) for
discussion on how to reduce these for particular boundary functions for reduced computational cost.
"""
function default_fixed_layers(
    constraint::UpperOnly,
    FT::Type{<:AbstractFloat},
)::Chain
    get_relus = Matrix{FT}([1 0; -1 0; 1 -1])   # (relu(f), relu(-f), relu(f - p))
    final_mul = Matrix{FT}([1 -1 -1])           # outputs relu(f) - relu(-f) - relu(f-p) = f - relu(f-p) = min(f, p)
    return Chain(
        apply_relus = Dense(get_relus, false, relu),
        sumlayer = MulLayer(final_mul),
    )
end

"""
    default_fixed_layers(constraint::TwoSided, FT::Type{<:AbstractFloat})::Chain

Creates the default fixed-weight layers for applying a lower and upper bound to an output. The general matrices for
arbitrary boundaries are used on output [f_upper, f_lower, p]ᵀ, where f_upper is the upper bound function output,
f_lower is the lower bound function output, and p is the output of the predictive model. See the NeuralDepthModel
paper (https://doi.org/10.1175/AIES-D-24-0040.1) for discussion on how to reduce these for particular boundary
functions for reduced computational cost.
"""
function default_fixed_layers(
    constraint::TwoSided,
    FT::Type{<:AbstractFloat},
)::Chain
    scale_and_relus = Matrix{FT}([1 0 0; -1 0 0; 0 1 0; 0 -1 0; 1 0 -1])
    get_min = Matrix{FT}([1 -1 -1 1 -1; 0 0 1 0 0; 0 0 0 1 0])
    final_mul = Matrix{FT}([1 1 -1])
    return Chain(
        apply_relus = Dense(scale_and_relus, false, relu),
        apply_upper = Dense(get_min, false, relu),
        sumlayer = MulLayer(final_mul),
    )
end

"""
   ConstrainedNeuralModel(
        ::Type{FT},
        pred;
        upper_bound = nothing,
        lower_bound = nothing,
        in_scales::Union{Nothing, Vector{FT}} = nothing,
        out_scale::Union{Nothing, FT} = nothing,
        fixed_layers::Union{Nothing, Chain} = nothing,
        trainable_constraints::Bool = false
    ) where {FT <: AbstractFloat} 

Constructs a ConstrainedNeuralModel type, from a provided predictive model and any specified bounds. The options to supply a
vector of input scaling factors, a known model output scaling, or custom fixed layers are provided, as well as indicate whether
the parameters of the supplied constraints are trainable or not. If no fixed layers are supplied, the appropriate default 
fixed layers for the supplied bounds are used. See the NeuralDepthModel paper (https://doi.org/10.1175/AIES-D-24-0040.1) 
for discussion on how to reduce these for particular boundary functions for reduced computational cost.

The predictive model need not be a Flux Layer/Chain or even a neural network, though Flux contains many optimizations for models defined as such.
Just ensure the passed predictive model is of a callable type (a functor), with a defined method for (model::YourModelType)(x::YourInput)
as well as a return for Flux.trainable(model::YourModelType). This enables the calibration of many types of parameterized models
while obeying prescribed functional constraints. The ConstrainedNeuralModel type was designed for 1D neural networks taking vector
inputs, but can readily work for multidimensional inputs and outputs - just make sure to pass a predictive model which flattens the output
to a 1xN Matrix at the end (requiring a subsequent reshaping layer after the ConstrainedNeuralModel output), and that specified constraint
functions similarly map tensor or multidimensional inputs to a correspondingly-ordered 1xN Matrix. The ConstrainedNeuralModel can also act
as a layer within a larger model chain. Constraint functions can also be functions of additional data or parameters beyond those of the inputs
and predictive model outputs, as long as the passed constraint type is wrapped into a callable functor of such entities, instead of a
function taking more than two arguments (the predictive model output and the input) - if you also want the constraints to be trainable,
make sure the right returns are defined for Flux.trainable() of your constraint functors.

If supplying your own fixed layers, note the ordering of outputs for your constraint type - [f, p]ᵀ for single upper- or lower- bounds,
where f is the bound function output and p is the output of the predictive model, or [f_upper, f_lower, p]ᵀ for both bounds, where f_upper
and f_lower specify the upper and lower bounds, respectively. As the output-scaling is already internally handled, do not specify
custom fixed layers which already include output scaling; either specify the output scaling at construction or make use of the
`set_predictive_model_out_scale!()` function.
"""
function ConstrainedNeuralModel(
    ::Type{FT},
    pred;
    upper_bound = nothing, #type is arbitrary instead of ::Function, as GPU-performant functor types could be passed here instead
    lower_bound = nothing,
    in_scales::Union{Nothing, Vector{FT}} = nothing,
    out_scale::Union{Nothing, FT} = nothing,
    fixed_layers::Union{Nothing, Chain} = nothing,
    trainable_constraints::Bool = false,
) where {FT <: AbstractFloat}

    if !isnothing(fixed_layers)
        @assert all(
            l -> (l isa MulLayer) || (l isa Dense && l.σ === relu),
            fixed_layers.layers,
        ) """
        Supplied fixed layers must either consist entirely of Flux.Dense Layers
        with relu() activation, or ConstrainedNeuralModels.MulLayer layers.
        """
    end

    if trainable_constraints &&
       (upper_bound isa Function || lower_bound isa Function)
        error(
            """
      Constraints cannot be trainable if utilized bounds are functions. Please specify
      a valid functor F for bounds with declared Flux.trainable(::typeof(F)) to create
      trainable constraints.
      """,
        )
    end

    constraints = buildConstraints(upper_bound, lower_bound)
    #turn-off the bound compliance checking, leaving it optional to user:
    #= if any(has_evaluation_mode.(Ref(constraints), [:static, :dynamic]))
        below code...
    else
        error(
            """
      Provided constraint system is not sufficient for either static or dynamic usage.
      Make sure provided bounds permit consistent evaluation type. See tutorial on
      designing bounds for more details.
      """,
        )
    end
    =#
    use_scaling =
        isnothing(in_scales) ? NoScaling{FT}() : ConstScaling(FT.(in_scales))
    use_out =
        isnothing(out_scale) ? ScaleOutput(FT(1)) : ScaleOutput(FT(out_scale))
    use_layers =
        isnothing(fixed_layers) ? default_fixed_layers(constraints, FT) :
        convert_model(fixed_layers, FT)

    use_pred = convert_model(pred, FT)

    model = ConstrainedNeuralModel(
        use_pred,
        constraints,
        use_scaling,
        use_out,
        use_layers,
        deepcopy(use_layers[1].weight), #base fixed layer matrix for resets after training-based scaling
        isnothing(fixed_layers),
        Val{trainable_constraints}(),
    )
    return model
end

"""
    set_predictive_model_out_scale!(model::ConstrainedNeuralModel, scale::Real = model.out_scale[])

Sets the ouptut scaling constant (field `out_scale`) of a ConstrainedNeuralModel.
Such scaling can be used in scaled training for better convergence and performance (see `scale_model!()`).
Note: This method cannot be used for a fixed (static) model where weights/parameter values are final.
"""
function set_predictive_model_out_scale!(
    model::ConstrainedNeuralModel,
    scale::Real,
)
    FT = eltype(model.out_scale.sc)
    @assert scale != 0 """
    set_predictive_model_out_scale!(): Cannot set the predictive model's scaling constant to 0.
    """
    if typeof(model.fixed_layers[1].weight) <: SArray
        throw(
            ArgumentError(
                "Cannot set the output scaling of a fixed (static) model.",
            ),
        )
    end
    model.out_scale.sc[1] = FT(scale)
    return nothing
end

"""
    scale_model!(model::ConstrainedNeuralModel, mode::Symbol)

Handles scaling transformations of a ConstrainedNeuralModel to enable scaled training,
using the output scaling value specified (or set, see `set_predictive_model_out_scale!()`) for
the model. When using mode `:scaled_train` for a model with output scaling C, a constrained model
(where predictive model output is usually multiplied by C) instead multiplies the output by 1
and scales the constraints by 1/C, which can map outputs into a more conducive training range (the
utilized training data should then also have the targets multiplied by 1/C; see the tutorial).
After training, this scaling can be reset to the correct system-applicable scaling using this function
with mode `:reset`.
Note: This method cannot be used for a fixed (static) model where weights/parameter values are final.
We recommend saving a model (see: `save_model()`) after applying a :reset instead of saving the
scaled form.
"""
function scale_model!(model::ConstrainedNeuralModel, mode::Symbol)
    sc = model.out_scale.sc[1]
    if typeof(model.fixed_layers[1].weight) <: SArray
        throw(ArgumentError("Cannot scale a fixed (static) model."))
    end
    if mode == :scaled_train
        model.fixed_layers[1].weight .= model.initial_fixed_layer .* 1 / sc
    elseif mode == :reset
        model.fixed_layers[1].weight .= model.initial_fixed_layer
    else
        throw(
            ArgumentError(
                "scale_model!(): Possible modes are :reset or :scaled_train, mode `$(mode)` not recognized.",
            ),
        )
    end
    return nothing
end

"""
    build_model_API(m::ConstrainedNeuralModel)

Allows user to build the API of methods and types associated
with a ConstrainedNeuralModel that are not immediately available within
ClimaLand or this module. Could be used to pass as a custom metadata string
within `save_model()`.

Note: for this to work, one must make sure all method definitions of
associated bounds are defined with the `@bound` macro, and all affiliated
custom bound functor types are specified with the `@bound_type` macro.
"""
function build_model_API(m::ConstrainedNeuralModel)
    #first check if bounds are valid for building an API:
    bounds = get_bounds(m.constraints)
    if all(is_valid_bound, bounds)
        transfer_data = assess_model_transferability(m; get_api_data = true)
        api = build_API(transfer_data)
        return api
    else
        error(
            """
            This model's constraints must be specified with the @bound
            and/or @bound_type macro to use this functionality. Specify
            all bound methods with @bound and all custom types with @bound_type.
            """,
        )
    end
end

"""
    build_model_bound_documentation(m::ConstrainedNeuralModel)

Allows user to build the documentation regarding the bounds of a 
ConstrainedNeuralModel, to allow reproducibility if working with custom types
that are not already defined within ClimaLand or this module. Could be
used to pass as a custom metadata string within `save_model()`.

Note: for this to work, one must make sure all method definitions of
associated bounds are defined with the @bound macro, and all affiliated
custom bound functor types are specified with the @bound_type macro.
"""
function build_model_bound_documentation(m::ConstrainedNeuralModel)
    #first check if bounds are valid for building docs:
    bounds = get_bounds(m.constraints)
    if all(is_valid_bound, bounds)
        bound_docs = build_bound_docs(m)
        return bound_docs
    else
        error(
            """
            This model's constraints must be specified with the @bound
            and/or @bound_type macro to use this functionality. Specify
            all bound methods with @bound and all custom types with @bound_type.
            """,
        )
    end
end

"""
    save_model(
        model::ConstrainedNeuralModel,
        param_filepath::String
         model_struct_filepath::String;
         user_param_metadata::String = "",
         user_model_metadata::String = ""
    )

Saves the ConstrainedNeuralModel structure and parameters (separately) into JLD2 files, allowing models
to be saved and sent, or loaded into new systems. The splitting of parameters and structure into separate
files allows the reconstruction of the same model using slightly different parameters (for example, in the
case of fine-tuning a trained ConstrainedNeuralModel online within a larger climate simulation, where all
sub-model parameter updates are calculated togther). The saved metadata also specifies any modules beyond
ClimaLand and this module that are needed for the model to work.

Fixed (SArray parameters) models can be saved as well as more generic models, though we recommend
saving the more generic/dynamic (Array) form of a model and calling `make_static_model()` on the loaded
generic model, rather than saving the fixed model and calling `make_dynamic_model()` on the loaded form.

NOTE: ConstrainedNeuralModels and any relevant user-specified functions/types for the model's constraints
must be already defined in the importing codespace for a given model to load correctly. To aid with this,
one can use the `build_model_bound_documentation()` function to create metadata explicitly articulating 
the code of their constraints (marked via the @bound and/or @bound_type macros), so such code (as text)
could be pasted into the codepsace and permit model construction even if the original script that built
it has been lost. One can also use the `build_model_API()` functionality to identify any nested/custom
function calls or types (and their documentation) the model requires which do not exist in base julia,
this module, and the ClimaLand modules, though this will not include code syntax. However, since the API
stores documentation, we recommend circumventing this by including the code syntax or informative descriptions
in their documentation. The API, bound documentation, and any custom user metadata can be combined to supply
the string of either the `user_param_metadata` or `user_model_metadata` arguments.
"""
function save_model(
    model::ConstrainedNeuralModel,
    param_filepath::String,
    model_struct_filepath::String;
    user_param_metadata::String = "",
    user_model_metadata::String = "",
)
    ps, build_func = Flux.Optimisers.destructure(model)
    ps_metadata = build_parameter_metadata(
        build_func,
        user_param_metadata,
        model_struct_filepath,
    )
    jldsave(param_filepath, trainable_params = ps, metadata = ps_metadata)

    model_metadata =
        build_model_metadata(model, user_model_metadata, param_filepath)
    jldsave(
        model_struct_filepath,
        build_func = build_func,
        metadata = model_metadata,
        accept_num_params = build_func.length,
    )
    return nothing
end

"""
    load_model(filepath::String)

Loads the JLD2 file pertaining to a given ConstrainedNeuralModel structural or parameter file,
allowing one to inspect the file's metadata, or access the model's parameters or build-function
(which rebuilds the model once called with a supplied/loaded vector of model parameters).

NOTE: ConstrainedNeuralModels and any relevant user-specified functions/types for the model's constraints
must be already defined in the importing codespace for a given model to load correctly. To aid with this,
one can use the `build_model_bound_documentation()` function to create metadata explicitly articulating 
the code of their constraints (marked via the @bound and/or @bound_type macros), so such code (as text)
could be pasted into the codepsace and permit model construction even if the original script that built
it has been lost. One can also use the `build_model_API()` functionality to identify any nested/custom
function calls or types (and their documentation) the model requires which do not exist in base julia,
this module, and the ClimaLand modules, though this will not include code syntax. However, since the API
stores documentation, we recommend circumventing this by including the code syntax or informative descriptions
in their documentation. The API, bound documentation, and any custom user metadata can be combined to supply
the string of either the `user_param_metadata` or `user_model_metadata` arguments.
"""
function load_model(filepath::String)
    return JLD2.load(filepath)
end

"""
    load_model(strucure_filepath::String, param_filepath::String)

Returns the ConstrainedNeuralModel created from a structural file's build-function as well as a parameter
file's parameter vector. Dispatched form if two filepaths are supplied.
"""
function load_model(param_filepath::String, strucure_filepath::String)
    build_data = JLD2.load(strucure_filepath)
    param_data = JLD2.load(param_filepath)
    if haskey(param_data, "trainable_params")
        return load_model(param_data["trainable_params"], build_data)
    else
        error("File type not supported.\n")
    end
end

"""
    load_model(params::Vector{<:AbstractFloat}, data::Dict)

Returns the ConstrainedNeuralModel created from a structural file's build-function as well as a parameter
file's parameter vector. Dispatched form if a parameter vector and structural data object are supplied.
"""
function load_model(params::Vector{<:AbstractFloat}, data::Dict)
    if haskey(data, "build_func")
        build_func = data["build_func"]
    else
        error("Model data does not seem to be a supported format.")
    end
    @assert length(params) == data["accept_num_params"] """
    The passed parameters are not the size this model takes.
    This model takes $(data["accept_num_params"]) parameters,
    and $(length(params)) were provided.
    """
    try
        load_model(params, build_func)
    catch e
        print(e)
        error("""
    Inspect model metadata before attempting a new load:
    --------- model metadata: -----------------------------
    $(data["metadata"])
    """)
    end
end

"""
    load_model(params::Vector, build_func::Flux.Optimisers.Restructure)

Returns the ConstrainedNeuralModel created from a structural file's build-function as well as a parameter
file's parameter vector. Dispatched form if a parameter vector and explicit build-function are supplied.
"""
function load_model(params::Vector, build_func::Flux.Optimisers.Restructure)
    return build_func(params)
end

"""
    inspect_model_metadata(data::Dict)

Displays the metadata of a loaded ConstrainedNeuralModel structural file. 
"""
function inspect_model_metadata(data::Dict)
    if haskey(data, "metadata")
        print(data["metadata"])
    else
        error("Unsupported file type.")
    end
end

"""
    make_static_model(m::ConstrainedNeuralModel; skip_check = false)

Returns a static (fixed) form of a more generic ConstrainedNeuralModel, for more optimal
computational performance (using SArray parameters and inputs).

The default method mode requires the model's constraints to be capable for :static
evaluation mode, which are assessed prior to returning the static model. This requires
usage of the `@bound`/`@bound_type` macros around affiliated bounds and bound-types,
though this check can be skipped by specifying `skip_check = true`.
"""
function make_static_model(m::ConstrainedNeuralModel; skip_check = false)
    try
        if skip_check || has_evaluation_mode(m.constraints, :static)
            return Adapt.adapt_structure(SArray, m)
        end
    catch e
        error(
            """
  Uncertainty in constraints for making a static model for this model type.
  Make sure adequate methods for StaticArray SVector/SMatrix inputs have been defined
  with the `@bound` macro for all model bounds, or specify `skip_check = true` to skip this check.
      """,
        )
    end
end

"""
    make_dynamic_model(m::ConstrainedNeuralModel; skip_check = false)

Returns a more generic form of a static ConstrainedNeuralModel, in order to permit
further prototyping on an otherwise static model (e.g., a saved/loaded fixed model,
though we recommend saving the generic model and calling make_static_model after loading).

The default method mode requires the model's constraints to be capable for :dynamic
evaluation mode, which are assessed prior to returning the dynamic model. This requires
usage of the `@bound`/`@bound_type` macros around affiliated bounds and bound-types,
though this check can be skipped by specifying `skip_check = true`.
"""
function make_dynamic_model(m::ConstrainedNeuralModel; skip_check = false)
    try
        if skip_check || has_evaluation_mode(m.constraints, :dynamic)
            return Adapt.adapt_structure(Array, m)
        end
    catch e
        error(
            """
  Uncertainty in constraints for making a dynamic model for this model type.
  Make sure adequate methods for <:Array or AbstractArray inputs have been defined
  with the `@bound` macro for all model bounds, or specify `skip_check = true` to skip this check.
      """,
        )
    end
end

"""
    trainmodel!(
        model,
        x_train,
        y_train,
        loss;
        nepochs::Int = 100,
        nbatch = 64,
        opt = Flux.Optimisers.RMSProp(),
        verbose = false,
    )

An optional helper function for training a model, abstracting/wrapping away some of the training set-up process.
Requires any differentiable loss function defined with three arguments (model, input, target), any Flux-enabled
model, as well as iterable instance of input features as well as their corresponding targets (able to be passed
to Flux.DataLoader() as (x_train, y_train)). Optional arguments exist for specifying the number of epochs, batch
size, optmizer algorithm, as well as to print training loss during the training process.
"""
function trainmodel!(
    model,
    x_train,
    y_train,
    loss;
    nepochs::Int = 100,
    nbatch = 64,
    opt = Flux.Optimisers.RMSProp(),
    verbose = false,
)
    train_loader = Flux.DataLoader(
        (x_train, y_train),
        batchsize = nbatch,
        partial = false,
        shuffle = true,
    )
    opt_state = Flux.setup(opt, model)
    for epoch in 1:nepochs
        for (x, y) in train_loader
            Flux.train!(loss, model, [(x, y)], opt_state)
        end
        if verbose && (epoch % 10 == 0)
            print(
                "Epoch: ",
                epoch,
                " | training loss: ",
                loss(model, x_train, y_train),
                "\n",
            )
        end
    end
end

include("./model_utilities.jl")

end
