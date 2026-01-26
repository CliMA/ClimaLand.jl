module FastPowerMooncakeExt

using FastPower, Mooncake
Mooncake.@mooncake_overlay FastPower.fastpower(x,y) = x^y

end