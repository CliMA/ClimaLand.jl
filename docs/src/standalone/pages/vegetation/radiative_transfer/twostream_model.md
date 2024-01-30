## The Two-Stream Scheme
In order to treat the effects of multiple scattering by cloud particles, aerosols and air molecules, the two-stream approximations are employed in most shortwave radiation (i.e., solar, 300-2500 nm) schemes presently used in LSMs for numerical weather prediction and climate modelling. In two-stream approximations, the radiation field is divided into the direct solar beam, plus the diffuse solar radiation (i.e., radiation scattered at least once), and in two directions, downward and upward fluxes. The angular distribution of scattered radiation is not computed in any further detail, which means they are considered to be isotropic (Raisaenen, 2002).

The two-stream approximation, or scheme has been used to deal with radiative transfer in the atmosphere for many years. The basic procedure in applying it to vegetation is to expand a complex function in the control equations into Legendre functions and then truncate them to the first order closure to get a simple solution (Dai, 2007). After reviewing several variants of the two-stream approximation model in the calculation of atmospheric radiation, Meador (1980) presented a unified form of the variants and introduced a new and improved method.

Dickinson 1983 introduced this new two-stream method to estimate radiative transfer in a vegetated canopy, and Sellers 1985 used the two-stream approximation to calculate values of hemispheric canopy reflectance in the visible or photosynthecially active radiation (PAR) and near-infrared (NIR) wavelength intervals. The two-stream approximation treatment has been widely used in land surface process models until nowadays.
The approximation assumes that diffuse radiative fluxes are isotropic in the upward and downward directions. Supposing that the upper and lower leaf optical properties are identical, the two-stream approximation used to model radiative transfer in plant canopies is given in the following form:

```math
-\overline{\mu}(dI^{\uparrow})/dL + [1 - (1 - \beta)\omega]I^{\uparrow} - \omega \beta I^{\downarrow} = \omega \overline{\mu} K \beta_0 \exp{(-KL)},\\
-\overline{\mu}(dI^{\downarrow})/dL + [1 - (1 - \beta)\omega]I^{\downarrow} - \omega \beta I^{\uparrow} = \omega \overline{\mu} K (1-\beta_0) \exp{(-KL)}
```

where I↑ and I↓ are the upward and downward diffuse radiative fluxes normalized by the incident flux respectively, μ is the cosine of the zenith angle of the incident beam, K is the optical depth of direct beam per unit leaf area and is equal to G(μ)/μ, G(μ) is the relative projected area of leaf elements in the direction cos−1μ, μ is the average inverse diffuse optical depth per unit leaf area
and is equal to

```math
\int_{0}^{1}[\mu^{\prime}/G(\mu^{\prime})]d\mu^{\prime}
```
μ′ is the direction of scattered flux, ω is the scattering coefficient and is equal to ρleaf +τleaf , and L is the cumulative LAI. β and β0 are upscattering parameters for the diffuse and direct beams respectively. (See Sellers 1985 for details)

These equations can be solved as an exact solution with appropriate boundary conditions. For direct incident radiation, the appropriate top boundary condition is I↓ = 0 for L = 0, and the bottom boundary condition is I↑ = ρs[I↓ + exp (−kLT )] for L = LT , where ρs is the soil reflectance and LT is the total LAI. The corresponding solution yielded is then:

```math
I^{\uparrow} = \frac{h_1\exp{(-KL)}}{\sigma} + h_2\exp{(-hL)} + h_3\exp{(hL)},\\
I^{\downarrow} = \frac{h_4\exp{(-KL)}}{\sigma} + h_5\exp{(-hL)} + h_6\exp{(hL)}
```

For diffuse radiation, the appropriate top boundary condition is I↓ = 1 for L
= 0, and the bottom boundary condition is I↑ = ρsI↓ for L = LT. Then, the corresponding solution is:

```math
I^{\uparrow} = h_7\exp{(-hL)} + h_8\exp{(hL)},\\
I^{\downarrow} = h_9\exp{(-hL)} + h_{10}\exp{(hL)}
```

where coefficients such as σ and h1 to h10 are given in Sellers 1985. Note that there is an error in the expression for h4 in the appendix of Sellers 1985. The correct expression may be found in Sellers 1996.

The model has the following parameters:

| Output | Symmbol | Unit | Range |
| :---         |     :---:      |    :---:      |     :---:   |
| Absorbed Photosynthetically Active Radiation  | APAR   | μmol m⁻² s⁻¹  | 0-1500 |
| Absorbed Near-Infrared Radiation              | ANIR   | μmol m⁻² s⁻¹  | 0-1500 |

| Drivers | Symbol | Unit | Range |
| :---         |     :---:      |    :---:      |     :---:   |
| Photosynthetically Active Radiation | PAR | μmol m⁻² s⁻¹  | 0--1500 |
| Leaf Area Index   | LAI   | m² m⁻² | 0--10 |

| Parameters | Symbol | Unit | Range |
| :---         |     :---:      |    :---:      |     :---:   |
| Canopy PAR Reflectance | $\alpha\_PAR\_{leaf}$  | -  | 0.0--1.0 |
| Canopy NIR Reflectance | $\alpha\_NIR\_{leaf}$  | -  | 0.0--1.0 |
| Canopy PAR Transmittance | $\tau\_PAR\_{leaf}$  | -  | 0.0--1.0 |
| Canopy NIR Transmittance | $\tau\_NIR\_{leaf}$  | -  | 0.0--1.0 |
| Canopy Emissivity        | $ϵ\_canopy$         | -  | 0.0--1.0 |
| Clumping index | $Ω$  | -  | 0.0--1.0 |
| Zenith angle | $θ_s$  | rad | 0--π |

| Constants | Symbol | Unit | Value |
| :---         |     :---:      |    :---:      |     :---:   |
| Leaf angle distribution | $l_d$ | - | 0.5 |
| Typical wavelength per photon PAR | $\lambda\_\gamma\_PAR$ | m | 5e-7
| Typical wavelength per photon NIR | $\lambda\_\gamma\_NIR$ | m | 1.65e-6
