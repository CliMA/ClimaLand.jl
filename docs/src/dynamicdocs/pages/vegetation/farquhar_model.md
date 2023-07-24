# Farquhar Model
This section breaks down the Farquar model that describes the biochemical process of photosynthesis in plants as environmental conditions change.

The biochemical processes within a leaf determine the rate of photosynthesis, particularly the diffusion of CO$_2$ into the leaf, the assimilation of CO$_2$ during photosynthesis, and the transpiration of water vapor. It takes into account factors such as light intensity, temperature, and CO$_2$ concentration to estimate the rate at which plants convert light energy into chemical energy through photosynthesis.

GPP is calculated based on the biochemistry of C3 and C4 photosynthesis
to determine potential (unstressed by water availability) leaf-level photosynthesis. This is calculated in terms of two potentially-limiting rates:

some text

```@raw html
<iframe src="https://clima.westus3.cloudapp.azure.com/jsserve/leaf_An"
   style="height:1400px;width:100%;">
</iframe>
```

some text

```@raw html
<iframe src="https://clima.westus3.cloudapp.azure.com/jsserve/leaf_An_ci"
   style="height:1400px;width:100%;">
</iframe>
```


## Rubisco limited rate

```math
\begin{equation}
a_1(T, c_a, VPD) =
\begin{cases}
      V_{cmax}(T)  \frac{(c_i(T, c_a, VPD) - \Gamma^*(T))}{(c_i(T, c_a, VPD) + K_c(T)*(1+o_i/K_o(T)))} & \text{for C3}\\
      V_{cmax}(T) & \text{for C4}
\end{cases}
\end{equation}
```

The dependence on the atmospheric CO$_2$ concentration $c_a$ (mol/mol) and vapor pressure deficit $VPD$ arise in the expression for $c_i$,
```math
\begin{align}
    c_i(T, c_a, VPD) = \max{(c_a(1-1/m(VPD)), \Gamma^*(T)}),
\end{align}
```
where and $m$ is the Medlyn factor (see Stomatal Conductance).

We also have
```math
    \Gamma^*(T) = \Gamma^*_{25}\exp\left(\Delta H_{\Gamma^*}\frac{T - T_o}{T_o R T}\right),
```

where $\Delta H_{\Gamma^*}$ is the activation energy per mol for $\Gamma^*$.

## Light limited rate

```math
\begin{equation}
a_2 =
\begin{cases}
      J(T, PAR) (c_i - \Gamma^*)/4(c_i + 2  \Gamma^*) & \text{for C3}\\
      J(T, PAR) & \text{for C4}
\end{cases}       
\end{equation}
```

where J is the rate of electron transport, which has units of mol photon per m$^2$ per s. It depends on $PAR$ via $APAR$, as described below, and on $T$ via the dependence on $J_{max}$.

J is given by the root of the equation
```math
\begin{align}
    \theta_j J^2 - (I + J_{max}) J + I J_{max} &= 0 \nonumber \\
    I &= \frac{\phi}{2} (APAR) \nonumber \\
    J_{max}(T) &= V_{cmax}(T)\times e \exp\left(\Delta H_{J_{max}}\frac{T - T_o}{T_o R T}\right),\nonumber \\
J(T, PAR) &= \frac{(I + J_{max} - \sqrt{(I + J_{max})^2 - 4\theta_j I \times J_{max}}}{2\theta_j},
\end{align}
```
where $\phi = 0.6$ and $\theta_j = 0.9$ are the quantum yield of photosystem II and a curvature function (Bonan's book), and $\Delta H_{J_{max}}$ is the energy of activation of $J_{max}$.

The total net carbon assimilation (A$_n$, mol CO$_2$ m$^{-2}$ s$^{-1}$) is given by the weighted sum of C3 and C4 net carbon assimilation fractions following:
```math
\begin{align}
A_n(T, PAR, VPD, c_a) = \text{max}(0, \text{min}(a_1 \beta, a_2) - R_d)
\end{align}
```

where $\beta$ is the moisture stress factor which is related to the mean soil moisture concentration in the root zone and R$_d$ is the leaf dark respiration calculated as 
```math
\begin{align}
    R_{d,25}(\psi_l) &= f V_{cmax,25}\beta(\psi_l), \nonumber \\
    R_d (T, \psi_l) & = R_{d,25}(\psi_l)\exp\left(\Delta H_{R_{d}}\frac{T - T_o}{T_o R T}\right),
\end{align}
```

where $f = 0.015$ is a constant, $\Delta H_{R_d}$ is the energy of activation for $R_d$, and finally 
Vcmax is calculated as 
```math
\begin{equation}
V_{cmax}(T) = V_{cmax,25} \exp\left(\Delta H_{Vcmax}\frac{T - T_o}{T_o R T}\right)\\
\end{equation}
```
with $V_{cmax,25}$ is a parameter (Vcmax at the reference temperature 25 C), and $\Delta H_{Vcmax} = 65,330 J/mol$.

The moisture stress factor is related to the leaf water potential $\psi_l$ as
```math
\begin{align}
    \beta = \frac{1+ \exp{(s_c \psi_c)}}{1+ \exp{(s_c(\psi_c - \psi_l))}},
\end{align}
```
where $s_c = 4$MPa$^{-1}$, $\psi_c = -2$MPa, and $\psi_l$ is the leaf water potential computed by the plant hydraulics model.

GPP is the total canopy photosynthesis calculated as the integral of leaf-level photosynthesis over the entire canopy leaf area index:
```math
\begin{align}
GPP(T, PAR, c_a, VPD, \theta_s) = A_n  (1 - \exp(-K LAI \Omega))/K.
\end{align}
```
This is not currently needed by other components, but is used for offline validation of the model.

We need to supply the following parameters and “drivers"

- ``K_{c,25}`` and $K_{o,25}$, $V_{cmax, 25}$, $\Gamma^*_{25},\phi$, $\theta_j$, $o_i$, $s_c$, $\psi_c$
- ``\psi_l``, to compute $\beta$
- Temperature $T$, $PAR$, $c_a$, VPD, $\theta_s$.
- 
| Output | Symbol | Unit | Range |
| :---         |     :---:      |    :---:      |     :---:   |
| Total net carbon assimilation | $A_n$   | mol CO$_2$ m$^{-2}$ s$^{-1}$  | 0--25 |

| Drivers | Symbol | Unit | Range |
| :---         |     :---:      |    :---:      |     :---:   |
| Photosynthetically Active Radiation | PAR | μmol m⁻² s⁻¹  | 0--1500 |
| Temperature | $T$  | °C  | 0--50 |

| Parameters | Symbol | Unit | Range |
| :---         |     :---:      |    :---:      |     :---:   |
| Moisture stress | $β$  | -  | 0-1 |
| Leaf Area Index   | LAI   | m² m⁻² | 1--10 |
| $CO_2$ concentration | $c_a$   | ppm | 300e-6--500e-6 |
| Vapor pressure deficit | VPD | Pa  | 500-10000 |
  
| Constants | Symbol | Unit | Value |
| :---         |     :---:      |    :---:      |     :---:   |
| Zenith angle | $θ_s$  | rad | 0.6 |
| Leaf angle distribution | $l_d$ | - | 0.5 |
| Canopy reflectance | $ρ_{leaf}$  | -  | 0.1 |
| Clumping index | $Ω$  | -  | 0.69 |
| $CO_2$ compensation at 25°C | Γ$^*_{25}$  | mol/mol | 4.275e-5 |
| A constant energy of activation | $ΔH_{Γ^*}$ | J/mol | 43540 |
| Standard temperature | $T_o$  | K | 298.15 |
| Universal gas constant | $R$  | J/mol | 8.314 |
| The maximum rate of carboxylation of Rubisco | $V_{cmax25}$  | mol CO$_2$ m$^{-2}$ s$^{-1}$ | 5e-5 |
| Energy of activation of $J_max$ | $ΔH_{J_max}$ | J/mol | 43540 |
| An empirical "curvature parameter" | $θ_j$  | -  | 0.9 |
| Clumping index | $Ω$  | -  | 0.69 |
| $CO_2$ compensation at 25°C | Γ$^*_{25}$  | mol/mol | 4.275e-5 |
| A constant energy of activation | $ΔH_{Γ^*}$ | J/mol | 43540 |
| Standard temperature | $T_o$  | K | 298.15 |
| Universal gas constant | $R$  | J/mol | 8.314 |
