module CompositionsBaseInverseFunctionsExt
using CompositionsBase
import InverseFunctions

InverseFunctions.inverse(::typeof(deopcompose)) = splat(opcompose)
InverseFunctions.inverse(::typeof(splat(opcompose))) = deopcompose
InverseFunctions.inverse(::typeof(decompose)) = splat(compose)
InverseFunctions.inverse(::typeof(splat(compose))) = decompose

end
