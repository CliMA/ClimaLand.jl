# PolygonOps.jl

```@index
```

# `inpolygon`

```@docs
PolygonOps.inpolygon
PolygonOps.HaoSun
PolygonOps.HormannAgathos
```

## Examples

```
using PolygonOps
using StaticArrays

xv = [ 0.05840, 0.48375, 0.69356, 1.47478, 1.32158, 
        1.94545, 2.16477, 1.87639, 1.18218, 0.27615, 
        0.05840 ]
yv = [ 0.60628, 0.04728, 0.50000, 0.50000, 0.02015, 
        0.18161, 0.78850, 1.13589, 1.33781, 1.04650, 
        0.60628 ]

xa = 0:0.1:2.3
ya = 0:0.1:1.4

polygon = SVector.(xv,yv)
points = vec(SVector.(xa',ya))

inside = [inpolygon(p, polygon; in=true, on=false, out=false) for p in points]

using Plots
plot(Tuple.(polygon), legend=false)
scatter!(Tuple.(points), marker_z=inside)
```


# `PolygonOps.area`

```@docs
PolygonOps.area
```
