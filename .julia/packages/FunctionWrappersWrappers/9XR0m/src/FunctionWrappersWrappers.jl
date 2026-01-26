module FunctionWrappersWrappers

using FunctionWrappers

export FunctionWrappersWrapper

struct FunctionWrappersWrapper{FW,FB}
  fw::FW
end
(fww::FunctionWrappersWrapper{FW,FB})(args::Vararg{Any,K}) where {FW,K,FB} = _call(fww.fw, args, fww)

_call(fw::Tuple{FunctionWrappers.FunctionWrapper{R,A},Vararg}, arg::A, fww::FunctionWrappersWrapper) where {R,A} = first(fw)(arg...)
_call(fw::Tuple{FunctionWrappers.FunctionWrapper{R,A1},Vararg}, arg::A2, fww::FunctionWrappersWrapper) where {R,A1,A2} = _call(Base.tail(fw), arg, fww)

const NO_FUNCTIONWRAPPER_FOUND_MESSAGE = "No matching function wrapper was found!"

struct NoFunctionWrapperFoundError <: Exception end

function Base.showerror(io::IO, e::NoFunctionWrapperFoundError)
    print(io, NO_FUNCTIONWRAPPER_FOUND_MESSAGE)
end

_call(::Tuple{}, arg, fww::FunctionWrappersWrapper{<:Any,false}) = throw(NoFunctionWrapperFoundError())
_call(::Tuple{}, arg, fww::FunctionWrappersWrapper{<:Any,true}) = first(fww.fw).obj[](arg...)

function FunctionWrappersWrapper(f::F, argtypes::Tuple{Vararg{Any,K}}, rettypes::Tuple{Vararg{DataType,K}}, fallback::Val{FB}=Val{false}()) where {F,K,FB}
  fwt = map(argtypes, rettypes) do A, R
    FunctionWrappers.FunctionWrapper{R,A}(f)
  end
  FunctionWrappersWrapper{typeof(fwt),FB}(fwt)
end

end
