# Running ClimaLand at FLUXNET sites

To run a FLUXNET site, we just need the site ID.
Currently, available sites are: Harvard ("US-Ha1"), Ozark ("US-MOz"), Niwot Ridge ("US-NR1"), Vaira Ranch ("US-Var"). 
More sites will be added soon.
For example to run the Ozark site:

```jl
julia> fluxnet_simulation("US-MOz")
```

This will run ClimaLand with default parameters and domain. Optionally, you can modify these:

```jl
julia> fluxnet_simulation("US-MOz"; params = ozark_default_params(; hetero_resp = hetero_resp_ozark(; b = 2)))
```

To get a list of ID of available FLUXNET sites, use:

```jl
julia> list_fluxnet_ID() # not implemented yet
```

To copy the parameters or domain template of a fluxnet site, use:

```jl
julia> template_parameters(site_ID) # not implemented yet
julia> template_domain(site_ID) # not implemented yet
```

# Calibrating parameters

The default parameters sites are already calibrated, but users can modify which parameter to calibrate, or 
during which period, priors, etc. (see full documentation). 

```jl
julia> calibrate(site_ID, parameters_to_calibrate) # not implemented yet
```
 
