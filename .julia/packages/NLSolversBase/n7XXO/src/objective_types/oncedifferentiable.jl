# Used for objectives and solvers where the gradient is available/exists
mutable struct OnceDifferentiable{TF, TDF, TX} <: AbstractObjective
    f # objective
    df # (partial) derivative of objective
    fdf # objective and (partial) derivative of objective
    F::TF # cache for f output
    DF::TDF # cache for df output
    x_f::TX # x used to evaluate f (stored in F)
    x_df::TX # x used to evaluate df (stored in DF)
    f_calls::Vector{Int}
    df_calls::Vector{Int}
end

### Only the objective
# Ambiguity
OnceDifferentiable(f, x::AbstractArray,
                   F::Real = real(zero(eltype(x))),
                   DF::AbstractArray = alloc_DF(x, F); inplace = true, autodiff = :finite,  
                   chunk::ForwardDiff.Chunk = ForwardDiff.Chunk(x)) =
    OnceDifferentiable(f, x, F, DF, autodiff, chunk)
#OnceDifferentiable(f, x::AbstractArray, F::AbstractArray; autodiff = :finite) =
#    OnceDifferentiable(f, x::AbstractArray, F::AbstractArray, alloc_DF(x, F))
function OnceDifferentiable(f, x::AbstractArray,
                   F::AbstractArray, DF::AbstractArray = alloc_DF(x, F);
                   inplace = true, autodiff = :finite)
    f! = f!_from_f(f, F, inplace)

    OnceDifferentiable(f!, x::AbstractArray, F::AbstractArray, DF, autodiff)
end


function OnceDifferentiable(f, x_seed::AbstractArray{T},
                            F::Real,
                            DF::AbstractArray,
                            autodiff, chunk) where T
    # When here, at the constructor with positional autodiff, it should already
    # be the case, that f is inplace.
    if  typeof(f) <: Union{InplaceObjective, NotInplaceObjective}

        fF = make_f(f, x_seed, F)
        dfF = make_df(f, x_seed, F)
        fdfF = make_fdf(f, x_seed, F)

        return OnceDifferentiable(fF, dfF, fdfF, x_seed, F, DF)
    else
        backend = get_adtype(autodiff, chunk)
        grad_prep = DI.prepare_gradient(f, backend, x_seed)
        function g!(_g, _x)
            DI.gradient!(f, _g, grad_prep, backend, _x)
            return nothing
        end
        function fg!(_g, _x)
            y, _ = DI.value_and_gradient!(f, _g, grad_prep, backend, _x)
            return y
        end
        return OnceDifferentiable(f, g!, fg!, x_seed, F, DF)
    end
end

has_not_dep_symbol_in_ad = Ref{Bool}(true)
OnceDifferentiable(f, x::AbstractArray, F::AbstractArray, autodiff::Symbol, chunk::ForwardDiff.Chunk = ForwardDiff.Chunk(x)) =
OnceDifferentiable(f, x, F, alloc_DF(x, F), autodiff, chunk)
function OnceDifferentiable(f, x::AbstractArray, F::AbstractArray,
                            autodiff::Bool, chunk::ForwardDiff.Chunk = ForwardDiff.Chunk(x))
    if autodiff == false
        throw(ErrorException("It is not possible to set the `autodiff` keyword to `false` when constructing a OnceDifferentiable instance from only one function. Pass in the (partial) derivative or specify a valid `autodiff` symbol."))
    elseif has_not_dep_symbol_in_ad[]
        @warn("Setting the `autodiff` keyword to `true` is deprecated. Please use a valid symbol instead.")
        has_not_dep_symbol_in_ad[] = false
    end
    OnceDifferentiable(f, x, F, alloc_DF(x, F), :forward, chunk)
end
function OnceDifferentiable(f, x_seed::AbstractArray, F::AbstractArray, DF::AbstractArray,
    autodiff::Symbol , chunk::ForwardDiff.Chunk = ForwardDiff.Chunk(x_seed))
    if  typeof(f) <: Union{InplaceObjective, NotInplaceObjective}
        fF = make_f(f, x_seed, F)
        dfF = make_df(f, x_seed, F)
        fdfF = make_fdf(f, x_seed, F)
        return OnceDifferentiable(fF, dfF, fdfF, x_seed, F, DF)
    else
        F2 = similar(F)
        backend = get_adtype(autodiff, chunk)
        jac_prep = DI.prepare_jacobian(f, F2, backend, x_seed)
        function j!(_j, _x)
            DI.jacobian!(f, F2, _j, jac_prep, backend, _x)
            return _j
        end
        function fj!(_y, _j, _x)
            y, _ = DI.value_and_jacobian!(f, _y, _j, jac_prep, backend, _x)
            return y
        end
        return OnceDifferentiable(f, j!, fj!, x_seed, F, DF)
    end
end

### Objective and derivative
function OnceDifferentiable(f, df,
                   x::AbstractArray,
                   F::Real = real(zero(eltype(x))),
                   DF::AbstractArray = alloc_DF(x, F);
                   inplace = true)


    df! = df!_from_df(df, F, inplace)

    fdf! = make_fdf(x, F, f, df!)

    OnceDifferentiable(f, df!, fdf!, x, F, DF)
end

function OnceDifferentiable(f, j,
                   x::AbstractArray,
                   F::AbstractArray,
                   J::AbstractArray = alloc_DF(x, F);
                   inplace = true)

    f! = f!_from_f(f, F, inplace)
    j! = df!_from_df(j, F, inplace)
    fj! = make_fdf(x, F, f!, j!)

    OnceDifferentiable(f!, j!, fj!, x, F, J)
end


### Objective, derivative and combination
function OnceDifferentiable(f, df, fdf,
    x::AbstractArray,
    F::Real = real(zero(eltype(x))),
    DF::AbstractArray = alloc_DF(x, F);
    inplace = true)

    # f is never "inplace" since F is scalar
    df! = df!_from_df(df, F, inplace)
    fdf! = fdf!_from_fdf(fdf, F, inplace)

    x_f, x_df = x_of_nans(x), x_of_nans(x)

    OnceDifferentiable{typeof(F),typeof(DF),typeof(x)}(f, df!, fdf!,
    copy(F), copy(DF),
    x_f, x_df,
    [0,], [0,])
end

function OnceDifferentiable(f, df, fdf,
                            x::AbstractArray,
                            F::AbstractArray,
                            DF::AbstractArray = alloc_DF(x, F);
                            inplace = true)

    f = f!_from_f(f, F, inplace)
    df! = df!_from_df(df, F, inplace)
    fdf! = fdf!_from_fdf(fdf, F, inplace)

    x_f, x_df = x_of_nans(x), x_of_nans(x)

    OnceDifferentiable(f, df!, fdf!, copy(F), copy(DF), x_f, x_df, [0,], [0,])
end
