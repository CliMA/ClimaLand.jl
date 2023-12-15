# Running ClimaLand at FLUXNET sites

To run a FLUXNET site, we just need the site ID, for example to run the Ozark site:

```jl
julia> fluxnet_simulation("US-MOz")
```

This will run ClimaLand with default parameters and domain. Optionally, you can modify these:

```jl
julia> fluxnet_simulation("US-MOz"; customparams = path_to_parameters, customdomain = path_to_domain)
```

where `customparams` and `customdomain` are modified from templates. Each site FLUXNET site has a default template. 

To get a list of ID of available FLUXNET sites, use:

```jl
julia> list_fluxnet_ID()
```

To copy the parameters or domain template of a fluxnet site, use:

```jl
julia> template_parameters(site_ID)
julia> template_domain(site_ID)
```

# Calibrating parameters

The default parameters sites are already calibrated, but users can modify which parameter to calibrate, or 
during which period, priors, etc. (see full documentation). 

```jl
julia> calibrate(site_ID, parameters_to_calibrate)
```
 
