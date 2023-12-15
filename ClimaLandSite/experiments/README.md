# Running ClimaLand at FLUXNET sites

To run a FLUXNET site, we just need the site ID, for example to run the Ozark site:

```jl
julia> run("US-MOz")
```

This will run ClimaLand with default parameters and domain. Optionally, you can modify these:

```jl
julia> run("US-MOz"; parameters = customparams, domain = customdomain)
```

where `customparams` and `customdomain` are modified from the templates provided in `path`

# Calibrating parameters

The default parameter sites are already calibrated, but users can modify which parameter to calibrate, or 
during which period, etc. (see full documentation). 

```jl
julia> calibrate("US-MOz", parameters_to_calibrate)
```
 
