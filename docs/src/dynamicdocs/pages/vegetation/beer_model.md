# Radiative transfer scheme
This section describes multiple models of radiative transfer 
through the vegetation canopy, implemented in ClimaLSM. 

## Beer's law
The absorbed PAR driving photosynthesis is calculated following the Beer-
Lambert law:

```math
APAR(PAR, \theta_s) = (PAR)(1 - \rho_{leaf})(1 - e^{(-K(\theta_s) LAI  \Omega)})
```

where PAR ≈ SW/2 is the incident moles of photons per meter squared per
second in the PAR window, approximated as half of the incident shortwave flux
if PAR is not directly available, ρleaf is the PAR canopy reflectance, K is the
vegetation extinction coefficient following Campbell (1998), LAI is the leaf area
index, θs is the zenith angle, and Ω is the clumping index following Braghiere
(2021). K, Ω and ρleaf are all unitless. LAI is in m² m⁻².
In order to compute K, we θs in radians and the leaf angle distribution ld
(unitless). K is then defined as

```math
K = l_d/\max{(\cos{(\theta_s)}, \epsilon)}
```

so that at night, when 3π/2 > θs > π/2, K is large (lots of extinction) and
non-negative. The small value ε prevents dividing by zero.

### Interactive APAR(PAR, LAI, ρleaf, k, Ω)

```@raw html
<iframe src="https://clima.westus3.cloudapp.azure.com/jsserve/beer_APAR"
   style="height:1400px;width:100%;">
</iframe>
```
