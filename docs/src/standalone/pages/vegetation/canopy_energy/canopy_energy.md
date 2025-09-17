# Canopy Energy Balance (Bigleaf)

This page documents the **canopy energy** module used in ClimaLand’s vegetation standalone runs.
The ClimaLand code can be found [here](https://github.com/CliMA/ClimaLand.jl/blob/main/src/standalone/Vegetation/canopy_energy.jl).

---

## Prognostic canopy temperature

The canopy temperature is modeled prognostically. Its time tendency is given by

```math
\frac{dT_c}{dt} =
\frac{
    LW_n + SW_n - (\mathrm{SHF} + \mathrm{LHF}) + F_\mathrm{roots}
}{
    a_c \cdot \max(A_\mathrm{leaf} + A_\mathrm{stem}, \epsilon)
},
```

where

 $T_{c}$ is the canopy temperature (K),

 $LW_{n}$ is the net long-wave radiation (W m$^{-2}$), positive if increasing canopy energy,

 $SW_{n}$ is the net short-wave radiation (W m$^{-2}$), positive if increasing canopy energy

 $\mathrm{SHF}$ is the sensible heat flux (W m$^{-2}$), positive if decreasing canopy energy,

 $\mathrm{LHF}$ is the latent heat flux (W m$^{-2}$), positive if decreasing canopy energy,

 $F_{\mathrm{roots}}$ is the energy flux associated with root–soil water transport (W m$^{-2}$), positive if increasing canopy energy,

 $a_{c}$ is the canopy heat capacity per unit area (J m$^{-2}$ K$^{-1}$),

 $A_{\mathrm{leaf}}$ and $A_{\mathrm{stem}}$ are the leaf and stem area indices (–),

 and $\epsilon$ is the machine epsilon for numerical stability.

---

## Root energy flux

- In **standalone mode** with prescribed ground conditions, this flux is set to zero.
- In coupled runs with a **prognostic soil model**, the flux is included to conserve energy, since the soil model accounts for the energy in soil water.

---
