# Radiative transfer scheme
This section describes multiple models of radiative transfer 
through the vegetation canopy, implemented in ClimaLand. 

## Beer's law
Plants utilize Photosynthetically Active Radiation (PAR) for the process of photosynthesis, during which they convert light energy into chemical energy, fueling the synthesis of sugars and other organic compounds. PAR refers to the portion of the electromagnetic spectrum that is essential for photosynthesis in plants. PAR includes wavelengths ranging from approximately 400 to 700 nanometers and corresponds to the visible light spectrum. The unit used to measure PAR is called micromoles per square meter per second (μmol/m²/s), representing the number of photons within the PAR range that strike a square meter of a surface per second.

The portion of PAR that is actually absorbed by the vegetation canopy for photosynthesis is called Absorbed Photosynthetically Active Radiation (APAR). The APAR driving photosynthesis is calculated following the Beer-
Lambert law:

```math
APAR(PAR, \theta_s) = (PAR)(1 - \rho_{leaf})(1 - e^{(-K(\theta_s) LAI  \Omega)})
```

where PAR ≈ SW/2 is the incident moles of photons per meter squared per
second in the PAR window, approximated as half of the incident shortwave flux.
If PAR is not directly available, $ρ_{leaf}$ is the PAR canopy reflectance, K is the
vegetation extinction coefficient following Campbell (1998), LAI is the leaf area
index, $θ_s$ is the zenith angle, and $Ω$ is the clumping index following Braghiere
(2021). $K$, $Ω$ and $ρ_{leaf}$ are all unitless. LAI is in m² m⁻².
In order to compute $K$, we need $θ_s$ in radians and the leaf angle distribution $l_d$
(unitless). K is then defined as

```math
K = l_d/\max{(\cos{(\theta_s)}, \epsilon)}
```

so that at night, when 3π/2 > $θ_s$ > π/2, $K$ is large (lots of extinction) and
non-negative. The small value ε prevents dividing by zero.

The model has the following parameters:

| Output | Symbol | Unit | Range |
| :---         |     :---:      |    :---:      |     :---:   |
| Absorbed Photosynthetically Active Radiation  | APAR   | μmol m⁻² s⁻¹  | 0-1500 |

| Drivers | Symbol | Unit | Range |
| :---         |     :---:      |    :---:      |     :---:   |
| Photosynthetically Active Radiation | PAR | μmol m⁻² s⁻¹  | 0--1500 |
| Leaf Area Index   | LAI   | m² m⁻² | 0--10 |

| Parameters | Symbol | Unit | Range |
| :---         |     :---:      |    :---:      |     :---:   |
| Canopy reflectance | $ρ_{leaf}$  | -  | 0.0--1.0 |
| Extinction coefficient  | $K$   | - | 0.0--1.0 |
| Clumping index | $Ω$  | -  | 0.0--1.0 |
| Zenith angle | $θ_s$  | rad | 0--π |
  
| Constants | Symbol | Unit | Value |
| :---         |     :---:      |    :---:      |     :---:   |
| Leaf angle distribution | $l_d$ | - | 0.5 |

### Interactive APAR(PAR, LAI, $ρ_{leaf}$, $K$, $Ω$)

```@raw html
<iframe src="https://clima.westus3.cloudapp.azure.com/jsserve/beer_APAR"
   style="height:1400px;width:100%;">
</iframe>
```
