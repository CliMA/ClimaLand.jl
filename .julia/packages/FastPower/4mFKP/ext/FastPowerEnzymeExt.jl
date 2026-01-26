module FastPowerEnzymeExt

using FastPower
import FastPower: fastpower
using Enzyme
using Enzyme.EnzymeRules: FwdConfig

function Enzyme.EnzymeRules.forward(config::FwdConfig,
        func::Const{typeof(FastPower.fastpower)},
        RT::Type{<:Union{Duplicated, DuplicatedNoNeed}},
        _x::Union{Const, Duplicated}, _y::Union{Const, Duplicated})
    x = _x.val
    y = _y.val
    ret = func.val(x, y)
    T = typeof(ret)
    if !(_x isa Const) 
        dxval = _x.dval * y * (fastpower(x,y - 1))  
    else 
        dxval = make_zero(_x.val)
    end
    if !(_y isa Const)  
        dyval = x isa Real && x<=0 ? Base.oftype(float(x), NaN) : _y.dval*(fastpower(x,y))*log(x)
    else 
        dyval = make_zero(_y.val)
    end 
    if RT <: DuplicatedNoNeed
        return convert(T,dxval + dyval)
    else
        return Duplicated(ret, convert(T, dxval + dyval))
    end
end

function EnzymeRules.augmented_primal(config::Enzyme.EnzymeRules.RevConfigWidth{1}, 
        func::Const{typeof(fastpower)}, ::Union{Type{<:Active}, Type{<:Const}}, x::Union{Const,Active}, y::Union{Const,Active})
    if EnzymeRules.needs_primal(config)
        primal = func.val(x.val, y.val)
    else
        primal = nothing
    end
    return EnzymeRules.AugmentedReturn(primal, nothing, nothing)
end

function EnzymeRules.reverse(config::Enzyme.EnzymeRules.RevConfigWidth{1}, 
        func::Const{typeof(fastpower)}, dret, tape, _x::Union{Const,Active}, _y::Union{Const,Active})
    x = _x.val
    y = _y.val
    dxval =  _x isa Const ? nothing : dret.val * y * (fastpower(x,y - 1))
    dyval = _y isa Const ? nothing : (x isa Real && x<=0 ? Base.oftype(float(x), NaN) : dret.val * (fastpower(x,y))*log(x))
    return (dxval, dyval)
end

end
