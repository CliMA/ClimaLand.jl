# ParamViz.jl - dynamic parameterisation web app 
![chrome_a0AHoCQMHV](https://github.com/CliMA/ParamViz.jl/assets/22160257/832adffe-5a5b-4d46-9d15-a088bcb4b460)
## Install ParamViz.jl: (unregistered for now)
```jl
julia> ]
pkg> add https://github.com/CliMA/ParamViz.jl
```
## Load packages:
```jl
julia> using ParamViz
julia> using Unitful: m, s, mol, μmol
julia> using ClimaLSM
julia> using ClimaLSM.Canopy
julia> FT = Float64
```
## Create structs:
```jl
julia> drivers = Drivers(("PAR (μmol m⁻² s⁻¹)", "LAI (m² m⁻²)"),
                         (FT.([0, 1500 * 1e-6]), FT.([0, 10])),
                         ((mol*m^-2*s^-1, μmol*m^-2*s^-1), (m^2*m^-2, m^2*m^-2))
                        )

julia> parameters = Parameters(("canopy reflectance, ρ_leaf",
                                "extinction coefficient, K",
                                "clumping index, Ω"),
                               (FT.([0, 1]), FT.([0, 1]), FT.([0, 1])),
                               ((m, m), (m, m), (m, m)) # dummy units, no conversion
                              )
julia> constants = Constants(("a", "b"), (FT(1), FT(2))) # dummy constants
julia> inputs = Inputs(drivers, parameters, constants)
julia> output = Output("APAR (μmol m⁻² s⁻¹)", [0, 1500 * 1e-6], (mol*m^-2*s^-1, μmol*m^-2*s^-1))
```
## Create method for parameterisation
```jl
julia> import ParamViz.parameterisation
julia> function parameterisation(PAR, LAI, ρ_leaf, K, Ω, a, b)   
         APAR = plant_absorbed_ppfd(PAR, ρ_leaf, K, LAI, Ω) 
         return APAR
       end
```
## Call webapp:
```jl
julia> webapp(parameterisation, inputs, output)
```
## Open app in your browser: 
Open your favorite browser and go to the URL http://localhost:9384/browser-display
