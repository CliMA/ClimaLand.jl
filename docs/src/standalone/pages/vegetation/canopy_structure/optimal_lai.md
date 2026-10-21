# Optimal LAI Model

The optimal LAI model, `ZhouOptimalLAIModel`, predicts the seasonal to decadal dynamics of the leaf area index (LAI) from optimality principles, following [Zhou2025](@citet). It combines the light use efficiency (LUE) model of photosynthesis of the P-model [Stocker2020](@cite) with an optimization of leaf area: the seasonal maximum LAI is the smaller of an energy-limited value, which maximizes the carbon gain of the canopy net of the cost of its leaves, and a water-limited value, at which the canopy transpires a fraction of the precipitation. Within the season, LAI follows the steady state of the daily potential GPP, with a lag. The model also sets the stem and root area indices and, through a competition between C3 and C4 plants, the C3 fraction of photosynthesis.

## Model Overview

At each timestep, the model:
1. computes the potential GPP of the canopy with the P-model, and updates the running totals of the potential GPP and of the climate (precipitation, potential evapotranspiration, vapor pressure deficit and air temperature)
2. sets the seasonal maximum LAI from its energy and water limits
3. computes the steady-state LAI from the daily potential GPP
4. relaxes LAI toward the steady-state LAI, which represents the lag of leaf growth
5. sets the leaf, stem and root area indices of the canopy, and its C3 fraction

LAI and the running totals are prognostic variables, advanced by the time-stepper; the other quantities are diagnosed from them at each timestep.

## Potential and Actual GPP

The model distinguishes two kinds of gross primary productivity (GPP):
- the potential GPP, $A_0$, which the canopy would reach if it absorbed all the incoming photosynthetically active radiation (PAR), that is, with a fraction of absorbed PAR (fAPAR) of 1, which would require an infinitely dense canopy
- the actual GPP of a canopy of finite LAI, $A = \text{fAPAR} \, A_0$

where fAPAR follows Beer's law, with the light extinction coefficient $k$:

```math
\text{fAPAR} = 1 - e^{-k \, \text{LAI}}
```

The potential GPP is the product of the photosynthetic photon flux density and of the light use efficiency of the P-model (Equations 1 and 2 of Zhou et al. 2025):

```math
A_0 = \text{LUE} \cdot \text{PPFD}, \qquad \text{LUE} = \phi_0 \, m_j \sqrt{1 - \left(\frac{c^*}{m_j}\right)^{2/3}}
```

where we have:
- the photosynthetic photon flux density at the top of the canopy, PPFD (mol photons m⁻² s⁻¹), from the downwelling PAR
- the intrinsic quantum yield, $\phi_0$ (mol CO₂ per mol photons), which depends on temperature
- the cost of the maximum rate of electron transport, $c^*$, with LUE = 0 where $m_j < c^*$
- the CO₂ limitation of electron transport, $m_j = (c_i - \Gamma^*)/(c_i + 2\Gamma^*)$ for C3 plants and 1 for C4 plants, with the CO₂ compensation point $\Gamma^*$ and the leaf-internal CO₂ partial pressure $c_i = \chi c_a$

The ratio of leaf-internal to ambient CO₂, $\chi$, is the optimum of the P-model at the air temperature and vapor pressure deficit (VPD) $D$:

```math
\chi = \frac{\xi + \sqrt{D} \, \Gamma^* / c_a}{\xi + \sqrt{D}}, \qquad \xi = \sqrt{\frac{\beta \, (K + \Gamma^*)}{1.6 \, \eta^*}}
```

with the ratio $\beta$ of the unit costs of carboxylation and transpiration, the Michaelis-Menten coefficient of Rubisco $K$, and the viscosity of water relative to 25 °C, $\eta^*$ (see the [P-model](../photosynthesis/pmodel.md) page). The potential GPP and $\chi$ of the canopy blend those of pure C3 and pure C4 canopies with the C3 fraction $f_{\mathrm{C3}}$ of the canopy (see C3/C4 Competition): $A_0 = f_{\mathrm{C3}} A_{0,\mathrm{C3}} + (1 - f_{\mathrm{C3}}) A_{0,\mathrm{C4}}$. The potential GPP uses the current VPD, while the $\chi$ of the seasonal maximum LAI (below) is evaluated at the growing-season VPD.

### Parameters of the Potential GPP

The potential GPP and $\chi$ use the unit cost ratios $\beta$ and the cost of electron transport $c^*$ of the optimal-LAI model (`optimal_lai_β_c3`, `optimal_lai_β_c4` and `optimal_lai_cstar`), not those of the photosynthesis model; the other P-model parameters, such as the quantum yield, are shared. The defaults, $\beta$ = 146 for C3 plants, 146/9 for C4 plants, and $c^*$ = 0.41, are those of pyrealm, with which Zhou et al. (2025) fitted their model and Lavergne et al. (2022) the C3/C4 competition. These fits, and the LAI parameters calibrated here (see Climate Tree Share), hold for the potential GPP and $\chi$ that these values give. The $\beta$ and $c^*$ of the photosynthesis model are instead calibrated against GPP and surface fluxes; their values shift the potential GPP and $\chi$, and so LAI and the C3/C4 balance, which have not been refitted to them. The two sets are not different plant traits: they keep the fits of the LAI model and of the competition valid when the GPP parameters are recalibrated.

### No Soil-Moisture Stress

The potential GPP carries no soil-moisture stress: its soil-moisture stress factor is 1. Zhou et al. (2025) applied the soil-moisture penalty of [Stocker2020](@citet), with the soil moisture of the SPLASH water balance [Davis2017](@cite), driven by the observed climate. Here the soil moisture is simulated, and depends on the transpiration of the simulated canopy, and so on LAI. A stressed potential GPP would make the seasonal maximum LAI depend on LAI itself, with a positive feedback where water limits LAI: a denser canopy dries the soil, which lowers the potential GPP and so raises the water-limited fAPAR, which is inversely proportional to the annual potential GPP (see Seasonal Maximum LAI). Without the stress, the supply of water limits LAI only through the precipitation that the canopy can transpire, $f_0 P_{\mathrm{annual}}$.

## Seasonal Maximum LAI

The seasonal maximum LAI follows from the smaller of an energy-limited and a water-limited fAPAR (Equations 11 and 12 of Zhou et al. 2025):

```math
\begin{aligned}
\text{fAPAR}_{\max} &= \min\left(\text{fAPAR}_{\mathrm{energy}},\ \text{fAPAR}_{\mathrm{water}}\right), \\
\text{fAPAR}_{\mathrm{energy}} &= 1 - \frac{z}{k \, A_{0,\mathrm{annual}}}, \\
\text{fAPAR}_{\mathrm{water}} &= \frac{c_a \, (1 - \chi)}{1.6 \, D_{\mathrm{gs}}} \, \frac{f_0 \, P_{\mathrm{annual}}}{A_{0,\mathrm{annual}}}, \\
\text{LAI}_{\max} &= -\frac{1}{k} \ln\left(1 - \text{fAPAR}_{\max}\right),
\end{aligned}
```

with $\text{fAPAR}_{\max}$ clamped to [0, 1), and where we have:
- the annual potential GPP, $A_{0,\mathrm{annual}}$ (mol CO₂ m⁻² yr⁻¹), the running yearly total of $A_0$ (see Prognostic Variables)
- the annual precipitation, rain and snow, $P_{\mathrm{annual}}$ (mol H₂O m⁻² yr⁻¹; 1 mm ≈ 55.5 mol H₂O m⁻²)
- the unit cost of building and maintaining leaves, $z$ (mol CO₂ m⁻² yr⁻¹), see Leaf Cost
- the ambient CO₂ partial pressure, $c_a$ (Pa), the CO₂ mixing ratio times the air pressure (400 ppm at 101325 Pa ≈ 40.5 Pa)
- the growing-season VPD, $D_{\mathrm{gs}}$ (Pa), see Growing-Season Climate
- the ratio of leaf-internal to ambient CO₂, $\chi$, at the VPD $D_{\mathrm{gs}}$: at the current air temperature in the steady-state LAI, and at the mean air temperature of the growing season, $T_{\mathrm{gs}}$, in the $\text{LAI}_{\max}$ that sets the tree retention and the stem and root area indices
- the fraction of the precipitation available for transpiration, $f_0$

The energy-limited fAPAR maximizes the net carbon profit of the canopy, $A_0 \, \text{fAPAR} - z \, \text{LAI}$ (Equations 7 to 10): leaves are added until the last one no longer gains more carbon than it costs, so a high leaf cost relative to the potential GPP reduces it. The water-limited fAPAR is the one at which the canopy transpires the fraction $f_0$ of the precipitation: transpiration, $1.6 \, G_s D$, and GPP, $G_s c_a (1 - \chi)$, share the stomatal conductance $G_s$ (Equations 4 and 5), so transpiring $f_0 P$ supports a GPP of $f_0 P \, c_a (1 - \chi) / (1.6 \, D)$.

The fraction of the precipitation available for transpiration depends on the aridity index, $\text{AI} = \text{PET}_{\mathrm{annual}} / P_{\mathrm{annual}}$ (Equation 6):

```math
f_0 = f_{0,\max} \exp\left[-0.604 \ln^2\left(\frac{\text{AI}}{1.9}\right)\right]
```

It peaks at $f_{0,\max}$ at the transition between energy- and water-limited climates (AI = 1.9) and falls toward both the arid and the humid extremes. As in Zhou et al. (2025), the potential evapotranspiration PET is that of Priestley and Taylor in SPLASH [Davis2017](@cite):

```math
\text{PET} = \frac{1.26}{\lambda_v M_w} \, \frac{\Delta}{\Delta + \gamma} \max(R_n, 0), \qquad R_n = (1 - 0.17) \, \text{SW}_d + \epsilon \left(\text{LW}_d - \sigma_{\mathrm{SB}} T^4\right)
```

where we have:
- the slope of the saturation vapor pressure curve, $\Delta$, and the psychrometric constant, $\gamma$ (Pa K⁻¹)
- the latent heat of vaporization, $\lambda_v$ (J kg⁻¹), and the molar mass of water, $M_w$ (kg mol⁻¹), which give PET in mol H₂O m⁻² s⁻¹
- the net radiation of the SPLASH reference surface, of albedo 0.17, $R_n$ (W m⁻²), from the downwelling shortwave and longwave radiation, $\text{SW}_d$ and $\text{LW}_d$, the canopy emissivity $\epsilon$, the Stefan-Boltzmann constant $\sigma_{\mathrm{SB}}$ and the air temperature $T$ (K)

As in SPLASH, negative (night-time) net radiation does not count, and the ground heat flux is neglected.

## Leaf Cost

The unit cost of leaves averages the costs of tree and grass leaves geometrically with the tree share of the vegetation, and scales its maintenance part with the season:

```math
z = z_{\mathrm{tree}}^{f_{\mathrm{tree}}} \, z_{\mathrm{grass}}^{1 - f_{\mathrm{tree}}} \left(1 - f_m + f_m \frac{M}{365\ \text{days}}\right)
```

where we have:
- the costs of tree and grass leaves for a year-round season at 25 °C, $z_{\mathrm{tree}}$ and $z_{\mathrm{grass}}$ (mol CO₂ m⁻² yr⁻¹); grasses have the higher cost, as $z$ includes the below-ground allocation that supplies the leaves
- the tree share of the vegetation, $f_{\mathrm{tree}}$, the rest, $1 - f_{\mathrm{tree}}$, being grasses
- the maintenance share of the cost, $f_m$, the rest, $1 - f_m$, being the cost of construction, which does not depend on the season
- the yearly maintenance load, $M$ (days), the running yearly total of the time above freezing, weighted by $Q_{10}^{(T - T_f - 25)/10}$ with the air and freezing temperatures $T$ and $T_f$ (K)

The maintenance load is 365 days for a year-round season at 25 °C, so the leaves of short, cold seasons are cheaper. The tree share is diagnosed from the simulated climate by default (`tree_share = PrognosticTreeShare()`, see Climate Tree Share), or prescribed as a scalar or a field, e.g. from the natural vegetation of the CLM surface data (`tree_share = clm_tree_share(surface_space)`).

## Steady-State LAI

The steady-state LAI, $L_s$, is the LAI in equilibrium with the GPP if the weather were held constant: $L_s = m A_s$, with the steady-state GPP $A_s = A_{0,\mathrm{daily}} (1 - e^{-k L_s})$, capped by $\text{LAI}_{\max}$. Its solution is (Equations 13 to 15 and 20 of Zhou et al. 2025):

```math
\begin{aligned}
L_s &= \min\left\{\mu + \frac{1}{k} W_0\left[-k\mu \, e^{-k\mu}\right],\ \text{LAI}_{\max}\right\}, \qquad \mu = m \, A_{0,\mathrm{daily}}, \\
m &= \frac{\sigma \, \text{GSL} \, \text{LAI}_{\max}}{A_{0,\mathrm{annual}} \, \text{fAPAR}_{\max}},
\end{aligned}
```

where we have:
- the daily potential GPP, $A_{0,\mathrm{daily}}$ (mol CO₂ m⁻² day⁻¹), the running daily total of $A_0$
- the principal branch of the Lambert W function, $W_0$, which satisfies $W(x) \, e^{W(x)} = x$
- the ratio of the steady-state LAI to the steady-state GPP, $m$
- the growing-season length, GSL (days), see Growing-Season Climate
- the departure, $\sigma$, of the seasonal cycle of LAI from a square wave that would hold $\text{LAI}_{\max}$ over the whole growing season

The steady-state LAI is zero unless $k\mu > 1$: it vanishes when the potential GPP is low, as in winter.

## LAI Dynamics

Zhou et al. (2025) represent the time lag of the allocation of carbon to leaves with a daily exponential moving average, $\text{LAI}_t = \alpha L_{s,t} + (1 - \alpha) \, \text{LAI}_{t-1}$ (Equation 16). Here LAI is a prognostic variable that relaxes toward $L_s$ continuously, with the same memory:

```math
\frac{d\,\text{LAI}}{dt} = \frac{L_s - \text{LAI}}{\tau_{\mathrm{LAI}}}, \qquad \tau_{\mathrm{LAI}} = \frac{1\ \text{day}}{\alpha}
```

The default $\alpha$ = 0.202 gives a response time of about 5 days; Zhou et al. (2025) use $\alpha$ = 0.067, about 15 days. The steady-state LAI is evaluated at every timestep, from the current running totals.

## Leaf Area Index of the Canopy

LAI follows the seasonal optimum, so it falls to zero with the steady-state LAI when the potential GPP is low. The leaf area index of the canopy, which radiative transfer, photosynthesis and the rest of the land model use, keeps for the tree share at least a fraction $f_r$ of $\text{LAI}_{\max}$ through the unfavourable season, the evergreen and semi-evergreen part of tree canopies:

```math
\text{LAI}_{\mathrm{canopy}} = \text{LAI} + f_{\mathrm{tree}} \max\left(f_r \, \text{LAI}_{\max} - \text{LAI},\ 0\right)
```

where $\text{LAI}_{\max}$ is evaluated with the $\chi$ of the mean air temperature of the growing season. As for prescribed LAI, values of $\text{LAI}_{\mathrm{canopy}}$ below 0.05 are set to zero. The evergreen share of trees is not modelled: it depends on biogeography and history more than on climate (spruce in boreal Canada, larch in eastern Siberia), and a share predicted from climate did not improve LAI over a constant one.

## Stem and Root Area Indices

By default, the stem and root area indices follow $\text{LAI}_{\max}$ (`SAI = AreaIndexFromMaxLAI{FT}(toml_dict, "SAI")`, and likewise for RAI), with factors per unit of $\text{LAI}_{\max}$ that the tree share weights between trees and grasses:

```math
\begin{aligned}
\text{SAI} &= \left[f_{\mathrm{tree}} \, s_{\mathrm{tree}} + (1 - f_{\mathrm{tree}}) \, s_{\mathrm{grass}}\right] \text{LAI}_{\max}, \\
\text{RAI} &= \left[f_{\mathrm{tree}} \, r_{\mathrm{tree}} + (1 - f_{\mathrm{tree}}) \, r_{\mathrm{grass}}\right] \text{LAI}_{\max}
\end{aligned}
```

Grasses have little stem area, as they have no wood, but more root area per unit of leaf area than trees. SAI is set to zero where the canopy is absent (LAI + SAI < 0.05). SAI and RAI can also be prescribed, as scalars or fields. The prescribed-LAI biomass model, given a maximum LAI such as that of MODIS, uses the same factors with the tree share of the CLM surface data. The canopy height and the rooting depth are those of the CLM surface data, as for prescribed LAI.

## Growing-Season Climate

The growing-season length, GSL (days), is the running yearly count of days above freezing. Zhou et al. (2025) use instead the length of the continuous period above 0 °C longer than 5 days.

The growing-season VPD, $D_{\mathrm{gs}}$, is the mean VPD of the moist growing season, when the air is above freezing and the water input of the last 30 days, $W_{30}$, is at least half the PET of the last 30 days, $\text{PET}_{30}$. Zhou et al. (2025) define the growing season by temperature alone; where it has a long dry season (seasonally dry tropics), the vegetation transpires mostly in the wet season, and the VPD of the dry season made the water limit too tight. Where the moist season is shorter than a month (deserts), $D_{\mathrm{gs}}$ blends into the mean VPD above freezing:

```math
D_{\mathrm{gs}} = w \, \overline{D}_{\mathrm{moist}} + (1 - w) \, \overline{D}_{\mathrm{warm}}, \qquad w = \min\left(\frac{n_{\mathrm{moist}}}{30\ \text{days}},\ 1\right)
```

where we have:
- the mean VPD of the moist growing season, $\overline{D}_{\mathrm{moist}}$, and above freezing, $\overline{D}_{\mathrm{warm}}$ (Pa)
- the yearly number of moist days, $n_{\mathrm{moist}}$ (days)

The water input $W$ is the precipitation $P$ that falls above freezing and the melt of a degree-day snow store $S$, which holds the precipitation that falls below freezing:

```math
\begin{aligned}
\frac{dS}{dt} &= \begin{cases} P, & T \le T_f, \\ -R_{\mathrm{melt}}, & T > T_f, \end{cases} \qquad
W = \begin{cases} 0, & T \le T_f, \\ P + R_{\mathrm{melt}}, & T > T_f, \end{cases} \\
R_{\mathrm{melt}} &= \min\left[F_{\mathrm{dd}} \, (T - T_f),\ S / (1\ \text{h})\right],
\end{aligned}
```

with the degree-day factor $F_{\mathrm{dd}}$ = 3 mm K⁻¹ day⁻¹. The climate tree share uses two more statistics of the growing season: its mean air temperature, $T_{\mathrm{gs}}$ (°C), the yearly sum of the degree-days above freezing divided by the yearly number of days above freezing; and its number of dry months, $n_{\mathrm{dry}}$, the yearly number of days above freezing that are not moist, in months (12/365 of the days).

## Climate Tree Share

By default (`tree_share = PrognosticTreeShare()`), the tree share is diagnosed from the simulated climate rather than prescribed by a map, so that it changes with the climate:

```math
f_{\mathrm{tree}} = \frac{1}{1 + \exp\left[-\left(b_0 + b_L \, L_{\mathrm{tree}} + b_d \, n_{\mathrm{dry}} + b_T \, T_{\mathrm{gs}}\right)\right]}
```

where we have:
- the seasonal maximum LAI that a C3 tree canopy would reach, $L_{\mathrm{tree}}$ (m² m⁻²): $\text{LAI}_{\max}$ with the annual potential GPP of C3 plants, the leaf cost of trees, $z_{\mathrm{tree}}$, scaled by the season (see Leaf Cost), and the $\chi$ of C3 plants at $T_{\mathrm{gs}}$ and $D_{\mathrm{gs}}$
- the number of dry months of the growing season, $n_{\mathrm{dry}}$, and its mean air temperature, $T_{\mathrm{gs}}$ (°C), see Growing-Season Climate

Trees need water through the dry season and enough productivity to pay for their canopy; grasses take over where the dry season is long. The coefficients, with the leaf costs, the maintenance share, $\sigma$ and the tree leaf retention, were calibrated in an offline emulator of the model (226 natural-vegetation points of a global 8° grid, 2008 ERA5) against MODIS LAI and the tree share of the natural vegetation in the CLM5 surface data. Nothing in $f_{\mathrm{tree}}$ depends on the simulated LAI, so there is no feedback through the leaf cost. In this mode, $f_{\mathrm{tree}}$ is also the tree share of the C3/C4 competition.

## C3/C4 Competition

The C3 fraction of photosynthesis comes from the competition between C3 and C4 plants of Lavergne et al. (2022), as implemented in pyrealm, rather than from the static map of the photosynthesis model. The proportional C4 advantage in annual potential GPP sets the expected C4 share of the open (non-tree) canopy, and C4 grasses are shaded out under trees:

```math
\begin{aligned}
f_{\mathrm{C4,open}} &= \frac{1}{1 + \exp\left[-k_{c34} \left(\Delta / e - q_{c34}\right)\right]}, \qquad \Delta = \frac{A_{0,\mathrm{C4,annual}} - A_{0,\mathrm{C3,annual}}}{A_{0,\mathrm{C3,annual}}}, \\
f_{\mathrm{C3g}} &= \left(1 - f_{\mathrm{C4,open}}\right) \left(1 - f_{\mathrm{tree}}\right), \qquad f_{\mathrm{C4g}} = f_{\mathrm{C4,open}} \left(1 - f_{\mathrm{tree}}\right), \\
f_{\mathrm{C3}} &= 1 - f_{\mathrm{C4g}},
\end{aligned}
```

where we have:
- the annual potential GPP of pure C3 and pure C4 canopies, $A_{0,\mathrm{C3,annual}}$ and $A_{0,\mathrm{C4,annual}}$ (mol CO₂ m⁻² yr⁻¹)
- the steepness and the midpoint of the logistic, $k_{c34}$ and $q_{c34}$
- the shares of C3 trees, C3 grasses and C4 grasses, $f_{\mathrm{tree}}$, $f_{\mathrm{C3g}}$ and $f_{\mathrm{C4g}}$, which sum to one
- the C3 fraction of photosynthesis, $f_{\mathrm{C3}}$, as trees are C3

They are shares of productivity, not of ground area. pyrealm divides the advantage by $\exp[1/(1 + \text{TC})]$, with TC an observed tree cover; with no such input, TC = 0 and the divisor is $e$.

With a prescribed tree share, which then sets only the leaf cost, the tree retention and the stem and root area indices, the tree share of the competition is the C3 tree cover of Lavergne et al. (2022), estimated from the realized C3 GPP:

```math
f_{\mathrm{tree}} = \min\left[\max\left(\frac{\mathrm{tc}(g)}{\mathrm{tc}(g_{\mathrm{ref}})},\ 0\right),\ 1\right], \qquad \mathrm{tc}(g) = a_{tc} \, g^{b_{tc}} + c_{tc}, \qquad g = \text{GPP}_{\mathrm{C3,annual}} \, M_C \, \frac{365\ \text{days}}{\text{GSL}}
```

where we have:
- the realized C3 GPP, $\text{GPP}_{\mathrm{C3,annual}}$ (mol CO₂ m⁻² yr⁻¹), the running yearly total of $A_{0,\mathrm{C3}} \, (1 - e^{-k \, \text{LAI}}) \, \beta_m$, with the soil-moisture stress factor of the canopy, $\beta_m$
- the molar mass of carbon, $M_C$ (kg mol⁻¹), which gives $g$ in kg C m⁻² yr⁻¹
- the coefficients of the tree-cover relation, $a_{tc}$, $b_{tc}$ and $c_{tc}$, and the C3 GPP at canopy closure, $g_{\mathrm{ref}}$

The relation was fitted where the growing season lasts all year; scaling the C3 GPP to a year-long season lets boreal forests, whose annual GPP is low only because their season is short, keep the tree cover of their productivity during it.

## Prognostic Variables

LAI and the running totals of the potential GPP and of the climate are prognostic variables, in `Y.canopy.biomass`, advanced by the time-stepper (time-integrated variables, `ClimaLand.TimeIntegratedVariable`). A running total $X$ of a rate $f$ over a window $\tau$ obeys

```math
\frac{dX}{dt} = \frac{f \, \tau - X}{\tau_s},
```

so that it tends to $\tau \bar{f}$, the total of $f$ over the window, with $\bar{f}$ the mean of $f$ over the smoothing timescale $\tau_s \ge \tau$. A smoothing timescale longer than a year filters the seasonal cycle out of the yearly totals. The sums marked "growing" have a memory that grows with the time since the start of the simulation, $\tau_s = \min(t_{\mathrm{age}} + 1\ \text{day},\ \tau_{\mathrm{long}})$, so that they average their whole history until it reaches $\tau_{\mathrm{long}}$ and hardly depend on their initial values.

| Variable | Rate $f$ | Window $\tau$ | Smoothing $\tau_s$ | Description |
| :--- | :--- | :---: | :---: | :--- |
| `LAI` | - | - | - | LAI, which relaxes toward $L_s$ over $\tau_{\mathrm{LAI}}$ |
| `A0_daily` | $A_0$ | 1 day | 3 days | Daily potential GPP, $A_{0,\mathrm{daily}}$ |
| `A0_annual` | $A_0$ | 1 year | $\tau_{\mathrm{long}}$ | Annual potential GPP, $A_{0,\mathrm{annual}}$ |
| `A0c3_annual`, `A0c4_annual` | $A_{0,\mathrm{C3}}$, $A_{0,\mathrm{C4}}$ | 1 year | $\tau_{\mathrm{long}}$ | Annual potential GPP of pure C3 and pure C4 canopies |
| `GPPc3_annual` | $A_{0,\mathrm{C3}} \, (1 - e^{-k \, \text{LAI}}) \, \beta_m$ | 1 year | $\tau_{\mathrm{long}}$ | Annual realized C3 GPP |
| `precip_annual` | $P$ | 1 year | $\tau_{\mathrm{long}}$ | Annual precipitation, $P_{\mathrm{annual}}$ |
| `PET_annual` | PET | 1 year | $\tau_{\mathrm{long}}$ | Annual potential evapotranspiration, $\text{PET}_{\mathrm{annual}}$ |
| `growing_days` | 1/day above freezing | 1 year | $\tau_{\mathrm{long}}$ | Growing-season length, GSL |
| `VPDgs_annual` | $D$ above freezing | 1 year | $\tau_{\mathrm{long}}$ | Sum of the VPD above freezing, which with GSL gives $\overline{D}_{\mathrm{warm}}$ |
| `water_30d`, `PET_30d` | $W$, PET | 30 days | 30 days | Water input and PET of the last 30 days, $W_{30}$ and $\text{PET}_{30}$ |
| `VPD_moist_annual`, `moist_days` | $D$ and 1/day, when moist | 1 year | growing | Sum of the VPD, and number of days, of the moist growing season |
| `degree_days`, `warm_days` | $T - T_f$ and 1/day, above freezing | 1 year | growing | Degree-days and number of days above freezing |
| `maintenance_days` | $Q_{10}^{(T - T_f - 25)/10}$/day, above freezing | 1 year | growing | Maintenance load, $M$ |
| `snow_store` | - | - | - | Snow store of the water input, $S$ |
| `age` | - | - | - | Time since the start of the simulation, $t_{\mathrm{age}}$ |

The diagnostics `a0d`, `a0a`, `a0c3`, `a0c4` and `pra` output `A0_daily`, `A0_annual`, `A0c3_annual`, `A0c4_annual` and `precip_annual` (the latter in m yr⁻¹).

### Initial Conditions

LAI and the running totals of the potential GPP, the precipitation, the PET, and the days and VPD above freezing start from the final state of a spin-up of the model (`experiments/long_runs/optimal_lai_spinup.jl`, at `ClimaLand.Artifacts.optimal_lai_state_path()`) or, where it has no data, from a climatology (`ClimaLand.Artifacts.optimal_lai_initial_conditions_path()`), with LAI from MODIS. The 30-day totals start at 30/365 of the yearly ones, and the snow store and the sums with a growing memory start at zero. See `ClimaLand.Simulations.set_canopy_component_initial_conditions!` for the details.

## Model Assumptions

1. **Water limitation enters once**: the potential GPP carries no soil-moisture stress, so the supply of water limits LAI only through the precipitation term of the water-limited fAPAR.
2. **Beer's law light extinction**: the absorption of light decays exponentially through the canopy.
3. **Optimal stomatal behavior**: the ratio of leaf-internal to ambient CO₂, $\chi$, is the optimum of the P-model.
4. **Climate inputs from the simulated climate**: the growing-season length and VPD, the aridity index behind $f_0$, the maintenance load and the climate tree share follow the simulated air temperature, humidity, radiation and precipitation, rather than maps.
5. **Continuous update**: LAI and the running totals are advanced every timestep by the time-stepper, not by a daily callback.

## Parameters

The default values are in `toml/default_parameters.toml`.

| Parameter | Symbol | Name in TOML | Unit | Default |
| :--- | :---: | :--- | :---: | :---: |
| Light extinction coefficient | $k$ | `optimal_lai_k` | - | 0.5 |
| Unit cost of tree leaves, for a year-round season at 25 °C | $z_{\mathrm{tree}}$ | `optimal_lai_z_tree` | mol CO₂ m⁻² yr⁻¹ | 9.92 |
| Unit cost of grass leaves, for a year-round season at 25 °C | $z_{\mathrm{grass}}$ | `optimal_lai_z_grass` | mol CO₂ m⁻² yr⁻¹ | 154 |
| Maintenance share of the leaf cost (the rest, $1 - f_m$, is construction) | $f_m$ | `optimal_lai_maintenance_share` | - | 0.43 |
| Temperature sensitivity of leaf maintenance | $Q_{10}$ | `optimal_lai_maintenance_q10` | - | 2 |
| Departure of the seasonal cycle of LAI from a square wave | $\sigma$ | `optimal_lai_sigma` | - | 1.09 |
| Smoothing factor, the inverse of the response time of LAI in days (about 5 days) | $\alpha$ | `optimal_lai_alpha` | - | 0.202 |
| Fraction of the precipitation available for transpiration at the energy–water transition | $f_{0,\max}$ | `optimal_lai_f0_max` | - | 0.65 |
| Fraction of $\text{LAI}_{\max}$ that trees keep through the unfavourable season | $f_r$ | `optimal_lai_tree_retention` | - | 0.5 |
| Smoothing timescale of the yearly running totals (2 years) | $\tau_{\mathrm{long}}$ | `optimal_lai_tau_long_term` | s | 6.3072 × 10⁷ |
| Unit cost ratios of the potential GPP, C3 and C4 | $\beta$ | `optimal_lai_β_c3`, `optimal_lai_β_c4` | - | 146, 16.2 |
| Cost of Jmax of the potential GPP | $c^*$ | `optimal_lai_cstar` | - | 0.41 |

The stem and root area indices per unit of $\text{LAI}_{\max}$, shared with the prescribed-LAI biomass model:

| Parameter | Symbol | Name in TOML | Unit | Default |
| :--- | :---: | :--- | :---: | :---: |
| Stem area index of trees per unit of $\text{LAI}_{\max}$ | $s_{\mathrm{tree}}$ | `SAI_coeff_tree` | - | 0.26 |
| Stem area index of grasses per unit of $\text{LAI}_{\max}$ | $s_{\mathrm{grass}}$ | `SAI_coeff_grass` | - | 0.05 |
| Root area index of trees per unit of $\text{LAI}_{\max}$ | $r_{\mathrm{tree}}$ | `RAI_coeff_tree` | - | 0.11 |
| Root area index of grasses per unit of $\text{LAI}_{\max}$ | $r_{\mathrm{grass}}$ | `RAI_coeff_grass` | - | 0.34 |

The climate tree share:

| Parameter | Symbol | Name in TOML | Unit | Default |
| :--- | :---: | :--- | :---: | :---: |
| Intercept of the logistic | $b_0$ | `optimal_lai_tree_b0` | - | -1.07 |
| Coefficient of the seasonal maximum LAI of a tree canopy | $b_L$ | `optimal_lai_tree_b_lai` | (m² m⁻²)⁻¹ | 0.0254 |
| Coefficient of the number of dry months of the growing season | $b_d$ | `optimal_lai_tree_b_dry` | month⁻¹ | -0.372 |
| Coefficient of the mean air temperature of the growing season | $b_T$ | `optimal_lai_tree_b_temp` | °C⁻¹ | 0.0917 |

The C3/C4 competition, with the coefficients fitted by Lavergne et al. (2022) and used by pyrealm:

| Parameter | Symbol | Name in TOML | Unit | Default |
| :--- | :---: | :--- | :---: | :---: |
| Steepness of the logistic of the C4 advantage | $k_{c34}$ | `optimal_lai_c3c4_k` | - | 6.63 |
| Midpoint of the logistic, the C4 advantage at which C3 and C4 are equally expected | $q_{c34}$ | `optimal_lai_c3c4_q` | - | 0.16 |
| Coefficient of the tree-cover relation | $a_{tc}$ | `optimal_lai_tc_a` | - | 15.60 |
| Exponent of the tree-cover relation | $b_{tc}$ | `optimal_lai_tc_b` | - | 1.41 |
| Offset of the tree-cover relation, so that the cover vanishes below a threshold GPP | $c_{tc}$ | `optimal_lai_tc_c` | - | -7.72 |
| C3 GPP at canopy closure | $g_{\mathrm{ref}}$ | `optimal_lai_tc_gpp_ref` | kg C m⁻² yr⁻¹ | 2.8 |

## Diagnostic Variables

The fields are in `p.canopy.biomass`.

| Variable | Symbol | Unit | Field | Diagnostic |
| :--- | :---: | :---: | :--- | :---: |
| Leaf area index of the canopy | $\text{LAI}_{\mathrm{canopy}}$ | m² m⁻² | `area_index.leaf` | `lai` |
| Stem area index | SAI | m² m⁻² | `area_index.stem` | `sai` |
| Root area index | RAI | m² m⁻² | `area_index.root` | - |
| Seasonal maximum LAI | $\text{LAI}_{\max}$ | m² m⁻² | `LAI_max` | - |
| Steady-state LAI | $L_s$ | m² m⁻² | `L_opt` | - |
| Potential GPP of the canopy, and of pure C3 and C4 canopies | $A_0$, $A_{0,\mathrm{C3}}$, $A_{0,\mathrm{C4}}$ | mol CO₂ m⁻² s⁻¹ | `OptVars.A0`, `OptVars.A0_c3`, `OptVars.A0_c4` | - |
| Ratio of leaf-internal to ambient CO₂ | $\chi$ | - | `OptVars.χ` | - |
| Fraction of the precipitation available for transpiration | $f_0$ | - | `f0` | - |
| Growing-season length | GSL | days | `GSL` | - |
| Growing-season VPD | $D_{\mathrm{gs}}$ | Pa | `vpd_gs` | - |
| Shares of C3 trees, C3 grasses and C4 grasses | $f_{\mathrm{tree}}$, $f_{\mathrm{C3g}}$, $f_{\mathrm{C4g}}$ | - | `composition.tree`, `composition.c3_grass`, `composition.c4_grass` | `ftr`, `fc3g`, `fc4g` |
| C3 fraction of photosynthesis | $f_{\mathrm{C3}}$ | - | 1 - `composition.c4_grass` | `fc3` |

With a prescribed tree share, `ftr` is the tree share of the C3/C4 competition, not the prescribed one.

## Implementation Notes

- The model is a biomass model, `ZhouOptimalLAIModel <: AbstractBiomassModel`: the leaf area index of the canopy is in `p.canopy.biomass.area_index.leaf`, as for `PrescribedBiomassModel`, and the biomass is updated first in `update_aux`, before radiative transfer and photosynthesis read the area indices and the C3 fraction.
- It requires the P-model of photosynthesis, `PModel`, whose parameters (with the $\beta$ and $c^*$ of the optimal-LAI model) and constants it uses for the potential GPP.
- The equations are in `src/standalone/Vegetation/optimal_lai.jl`, and the model, its cache and its tendencies in `src/standalone/Vegetation/biomass.jl`.
- Units: 1 mm of precipitation is 1 kg m⁻², or 55.5 mol H₂O m⁻²; 400 ppm of CO₂ at 101325 Pa is 40.5 Pa; the P-model gives the light use efficiency in kg C per mol photons, which the molar mass of carbon, 0.0120107 kg mol⁻¹, converts to mol CO₂ per mol photons.

## References

```@bibliography
Pages = [@__FILE__]
Canonical = false
```
