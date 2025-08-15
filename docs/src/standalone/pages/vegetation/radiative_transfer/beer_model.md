## Beer's law

Plants absorb, transmit, and reflect shortwave radiation; the fraction of downwelling
radiation partitioned into each of these
categories under Beer's law is given by per wavelength band as
```math
\rm{abs}_\lambda = (1 - \alpha_{\rm{leaf}, \lambda})(1 - e^{(-K(\theta_s) LAI  \Omega)})\\
\rm{trans}_\lambda = e^{(-K(\theta_s) LAI  \Omega)}\\
\rm{refl}_\lambda = 1 - \rm{abs}_\lambda - trans_\lambda * (1-\alpha_{\rm{ground}})
```

where λ reflects the wavelength of light, SW\_d is the downwelling radiative flux,
α\_{leaf,λ} is the albedo of the leaves in that wavelength band, K is the extinction
coefficient, θ\_s is the zenith angle, LAI is the leaf area index, Ω
is the clumping index, and α\_{ground,λ} is the ground albedo.

The extinction coefficient is defined as

```math
K = l_d/\max{(\cos{(\theta_s)}, \epsilon)}
```
where $l_d$ is the leaf distribution factor, and where the denominator is structured
so that at night, when 3π/2 > $θ_s$ > π/2, $K$ is large (lots of extinction) and
non-negative. The small value ε prevents dividing by zero.

The model has the following variables:

| Output | Symbol | Unit | Range |
| :---         |     :---:      |    :---:      |     :---:   |
| Absorbed fraction of radiative flux per band | abs_λ | W m⁻²  | 0--1 |
| Reflected fraction of radiative flux per band | refl_λ | W m⁻²  | 0--1 |
| Transmitted fraction of radiative flux per band | trans_λ | W m⁻²  | 0--1 |

| Input | Symbol | Unit | Range |
| :---         |     :---:      |    :---:      |     :---:   |
| Leaf reflectance (albedo) | $α\_{\rm{leaf},λ}$  | -  | 0--1 |
| Extinction coefficient  | $K$   | - | K>0 |
| Clumping index | $Ω$  | -  | 0--1 |
| Zenith angle | $θ_s$  | rad | 0--π |
| Leaf Area Index   | LAI   | m² m⁻² | 0--10 |
| Leaf angle distribution | $l_d$ | - | 0--1 |
| Ground albedo | $α\_{\rm{ground},λ}$  | -  | 0--1 |