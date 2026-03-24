## Beer's law

Plants absorb, transmit, and reflect shortwave radiation; On both the downwards and upwards
pass of radiation through the canopy, the canopy transmits a fraction per wavelength band of
```math
f_{\rm{trans}} = e^{(-K(\theta_s) LAI  \Omega)}.\\
```
The total fraction of the incident radiation on the land surface absorbed by the canopy is
```math
f_{\rm{abs, canopy}} = (1 - \alpha_{\rm{leaf}})(1 - f_{\rm{trans}})(1+\alpha_{\rm{ground}}f_{\rm{trans}_\lambda}),\\
```
while the total absorbed by the ground is
```math
f_{\rm{abs, ground}} = (1 - \alpha_{\rm{ground}})f_{\rm{trans}}.\\
```

The upwelling radiation from the land is
```math
f_{\rm{refl}} = 1 - f_{\rm{abs, canopy}} - f_{\rm{trans}} * (1-\alpha_{\rm{ground}}).
```

where α\_{leaf,λ} is the albedo of the leaves in that wavelength band, K is the extinction
coefficient, θ\_s is the zenith angle, LAI is the leaf area index, Ω
is the clumping index, and α\_{ground} is the ground albedo. Each of these is a function of wavelength.

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
| Absorbed fraction of radiative flux per band | f_abs_λ | W m⁻²  | 0--1 |
| Reflected fraction of radiative flux per band | f_refl_λ | W m⁻²  | 0--1 |
| Transmitted fraction of radiative flux per band | f_trans_λ | W m⁻²  | 0--1 |

| Input | Symbol | Unit | Range |
| :---         |     :---:      |    :---:      |     :---:   |
| Leaf reflectance (albedo) | $α\_{\rm{leaf},λ}$  | -  | 0--1 |
| Extinction coefficient  | $K$   | - | K>0 |
| Clumping index | $Ω$  | -  | 0--1 |
| Zenith angle | $θ_s$  | rad | 0--π |
| Leaf Area Index   | LAI   | m² m⁻² | 0--10 |
| Leaf angle distribution | $l_d$ | - | 0--1 |
| Ground albedo | $α\_{\rm{ground},λ}$  | -  | 0--1 |



Beer's law as implemented has the following simplifications:

1. No leaf transmittance: Assumes leaves only reflect ($α_{\rm{leaf}}$) or absorb. Does not account for light transmitted through leaves ($τ_{\rm{leaf}}$).
2. Single ground reflection: Accounts for one reflection from ground back through canopy, but not infinite cascade.
3. No within-canopy scattering: Light scattered by leaves is not tracked within the canopy volume.

For vegetation with significant leaf transmittance (especially NIR where $τ_{\rm{leaf}}$ ~ 0.4-0.5) or applications requiring accurate multiple scattering, use the TwoStreamModel instead.
