# Autotrophic Respiration

This page documents the autotrophic respiration model by [Clark2011](@citet), used in the JULES model.
The ClimaLand code can be found [here](https://github.com/CliMA/ClimaLand.jl/blob/main/src/standalone/Vegetation/autotrophic_respiration.jl#L240).

---

## Model formulation

Autotrophic respiration of the canopy is partitioned into maintenance and growth components.
At the canopy level, the respiration flux (mol CO$_2$ m$^{-2}$ s$^{-1}$) is

```math
R_a = \Bigl(R_{pm} + R_g\Bigr)\; \frac{1 - \exp\!\left(-K \, LAI \, \Omega\right)}{K \, \Omega},
```

where $R_{pm}$ is the maintenance respiration (mol CO$_2$ m$^{-2}$ s$^{-1}$), $R_g$ is the growth respiration (mol CO$_2$ m$^{-2}$ s$^{-1}$), $K$ is the canopy extinction coefficient (–), $LAI$ is the leaf area index (m$^2$ leaf m$^{-2}$ ground), and $\Omega$ is a clumping factor (–).

---

### Nitrogen allocation

Nitrogen is distributed among leaves, roots, and stems according to

```math
S_c = \eta_{sl}\, h\, LAI \, H(SAI), \qquad
R_c = \sigma_l \, RAI
```

```math
n_m = \frac{Vc_{max25}}{n_e}, \qquad
N_l = n_m \, \sigma_l \, LAI, \qquad
N_r = \mu_r \, n_m \, R_c, \qquad
N_s = \mu_s \, n_m \, S_c
```

where:

 $H$ is the Heaviside function (–),

 $Vc_{max25}$ is the leaf-level maximum carboxylation rate at 25 °C (mol CO$_2$ m$^{-2}$ s$^{-1}$),

 $n_e$ is the conversion factor from $Vc_{max25}$ to nitrogen (mol CO$_2$ m$^{-2}$ s$^{-1}$ kg C$^{-1}$),

 $\eta_{sl}$ is the live stem wood coefficient (kg C m$^{-3}$),

 $h$ is the canopy height (m),

 $\sigma_l$ is the specific leaf density (kg C m$^{-2}$ leaf),

 $\mu_r$ is the root-to-leaf nitrogen ratio (–),

 $\mu_s$ is the stem-to-leaf nitrogen ratio (–),

 $N_l, N_r, N_s$ are the nitrogen contents in leaves, roots, and stems (kg N m$^{-2}$ ground),

 $S_c$ is the stem carbon pool (kg C m$^{-2}$ ground),

 $R_c$ is the root carbon pool (kg C m$^{-2}$ ground),

 $n_m$ is the nitrogen per unit carboxylation capacity (kg N mol$^{-1}$ CO$_2$ m$^{2}$ s).

---

### Maintenance respiration

```math
R_{pm} = R_d \,\Bigl(\beta + \frac{N_r + N_s}{\max(N_l, \, \epsilon)}\Bigr),
```

where $R_d$ is a base respiration rate (mol CO$_2$ m$^{-2}$ s$^{-1}$), $\beta$ is a parameter scaling leaf respiration (–), and $\epsilon$ is machine epsilon (–) for numerical stability.

---

### Growth respiration

```math
R_g = R_{el} \, \Bigl(A_n - R_{pm}\Bigr),
```

where $R_{el}$ is the relative contribution of growth respiration (–), and $A_n$ is the net assimilation rate (mol CO$_2$ m$^{-2}$ s$^{-1}$).

---

## Notes

- This formulation follows the approach used in large-scale ecosystem models such as JULES (Clark et al., 2011).
- Growth respiration is proportional to net assimilation after subtracting maintenance costs.
- Maintenance respiration scales with tissue nitrogen pools.
