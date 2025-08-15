The P-model is model for photosynthesis and stomatal conductance that extends the Farquhar model in two main ways: 1) it assumes that all plants adjust their internal CO2 concentration such that carbon assimilation is maximized relative to the combination of two costs—the cost of supporting RuBisCo-limited photosynthesis and the cost of transpiration through stomata; 2) it assumes that RuBisCo-limited and light-limited assimilation rates are equal, i.e., that $A_c = A_j$ (coordination hypothesis). These two additional constraints allow for the prediction of photosynthetic parameters such as $V_{c\max}$, $J_{\max}$, and stomatal conductance $g_{s}$. 
## Usage in ClimaLand
The P-model differs from other canopy component models in two main ways: 
1) The P-model encompasses *both* a photosynthesis model *and* a stomatal conductance model. Therefore, `PModel` must be used for photosynthesis iff `PModelConductance` is used for stomatal conductance. 
2) The P-model requires an extra callback `PModelCallback` which checks if it is local noon every model timestep and updates the optimal parameters according to local noon (see Implementation for scientific background). When you are constructing a simulation via `SciMLBase`, this callback can be constructed like so:

```julia
# make the callback 
pmodel_cb = ClimaLand.make_PModel_callback(FT, start_date, t0, dt, canopy)

# add this callback to the CallbackSet with driver, diag
cb = SciMLBase.CallbackSet(driver_cb, diag_cb, pmodel_cb);
```

**To be implemented**: when using the P-model within a LandSimulation, automatically detect this and use the P-model callback. 
## Theory 
The first assumption, which we'll call the "cost minimization principle", can be mathematically written as follows. We define a cost which is dependent on the transpiration $E$, maximum rate of Rubisco limited photosynthesis $V_{c\max}$, normalized to the rate of assimilation $A$:
```math
$$ \mathrm{Cost} = a \frac{E}{A} + b \frac{V_{c\max}}{A} $$
```
$a$ and $b$ are dimensionless numbers that represent the relative weight of each cost. Minimization of this cost represents the general principle that plants seek to maximize assimilation while keeping water loss low ($E$ term) and not overly investing in synthesizing Rubisco protein. The parameter which is adjusted is $\chi = \dfrac{c_{i}}{c_{a}}$, the ratio of internal to ambient (external) CO2 concentration. Thus, the first constraint is mathematically written as
```math
$$ a \frac{\partial (E/A)}{\partial \chi} = -b \frac{\partial (V_{cmax} / A)}{ \partial \chi} $$
```
The second assumption, called the "coordination hypothesis" (see Chen et al., 1993), states that $V_{c\max}$
varies according to APAR such that it is neither in excess or in deficit of what is required for full utilization of the light. Mathematically, this means that $A_{c} = A_{j}$. These two constraints, when added to the Farquhar model, allows for the computation of $V_{c\max}$, $J_{\max}$, and $\xi$ (sensitivity parameter to VPD; see Wang 2017) from environmental conditions: 
```math
\begin{align*}
\xi^\mathrm{opt} &= \sqrt{ \dfrac{\beta(K + \Gamma^*)}{\eta^* D_{H_{2}O}/D_{CO_{2}}} } \\[0.2in]
J_{\max}^\mathrm{opt} &= \dfrac{4\phi_{0} \cdot  \mathrm{APAR}}{\sqrt{ \left[ \beta_{m}^2\left(1 - \left[ c^* \dfrac{c_{i} + 2\Gamma^*}{c_{i} - \Gamma^*} \right]^{2/3} \right) \right]^{-1} - 1 }} \\[0.2in]
V_{c\max}^{\mathrm{opt}} &= \phi_{0} \beta_{m} \cdot \mathrm{APAR} \cdot \dfrac{c_{i} + K}{c_{i} + 2\Gamma^*} \sqrt{ 1 - \left[ c^*  \dfrac{c_{i}+ 2\Gamma^*}{c_{i} - \Gamma^*}\right]^{2/3} }
\end{align*}
```
We refer the reader to the Farquhar docs page for definitions of the biochemical parameters referenced above.

## Implementation
We follow the scheme of Mengoli et al., (2022) who propose updating optimal parameters $V_{c\max 25}$, $J_{\max25}$, and $\xi$ according to conditions at local noon (near maximal APAR). These are then timestepped via the forward Euler discretization of the following ODE: 
```math
$$
\tau \dfrac{ d\bar{x}}{dt} = \bar{x} - x \implies \bar{x}_{t+1} = \alpha \bar{x}_{t} + (1 - \alpha) x_{t+1}
$$
```
where $\bar{x}$ is the acclimated parameter, $x$ is the parameter computed at local noon, $\alpha = 1 - \dfrac{\Delta t}{\tau}$ and $\tau$ is the timescale of acclimation. Since we update this equation every local noon, $\Delta t$ is 1 day and $\tau = 15$ days corresponds to $\alpha = \dfrac{14}{15} \approx 0.933$ (the default timescale used in Mengoli 2022). The optimal variables $\bar{x}$ are stored in the model cache in `p.canopy.photosynthesis.OptVars`.  

At every model timestep, the latest acclimated values $V_{c\max 25}^{\mathrm{opt}}$, $J_{\max 25}^{\mathrm{opt}}$, are then adjusted via modified-Arrhenius type functions of form 
```math
$$
V_{c\max} = V_{c\max 25}^{\mathrm{opt}} \cdot \underbrace{ \exp \left( \dfrac{\Delta H_{a}(T - T_{0})}{T_{0} TR} \right) }_{ \text{activation (Arrhenius)} } \cdot \underbrace{ \dfrac{1 + e^{(T_{0}\Delta S - \Delta H_{d})/RT_{0}}}{1 + e^{(T\Delta S - \Delta H_{d})/RT}}  }_{ \text{deactivation} }
$$
```
where $\Delta H_{a}$, $\Delta H_{d}$ are standard enthalpies of activation and deactivation, $T_{0} = 298.15$ K is the reference temperature, and $\Delta S$ is the entropy change in deactivation. The acclimated $\xi^\mathrm{opt}$ is used to compute the instantaneous intercellular CO2 concentration $c_{i}$
```math
$$
c_{i} =\dfrac{ \xi^\mathrm{opt} c_{a} + \Gamma^* \sqrt{ \max(D,0) }}{\xi^\mathrm{opt} + \sqrt{ \max(D,0) }}
$$
```
where $D$ is the vapor pressure deficit (VPD). Finally, carboxylation- and light-limited assimilation is computed as
```math
\begin{align*}
A_{c} &= V_{c\max} \dfrac{c_{i} - \Gamma^*}{c_{i} + K} \\
A_{j} &= \dfrac{1}{4} \underbrace{ \dfrac{4\phi_{0} I_{abs}}{\sqrt{ 1 + \left(\dfrac{4\phi_{0}I_{abs}}{J_{\max}}\right)^2 }} }_{ J } \cdot \dfrac{c_{i} - \Gamma^*}{c_{i} + 2 \Gamma^*}
\end{align*}
```
The remaining steps for computing $A_{n}$ and upscaling to canopy-level GPP are identical to the Farquhar model.

The P-model also gives a stomatal conductance model because it predicts $\chi = \dfrac{c_{a}}{c_{i}}$. From Fick's law, we have that 
```math
$$
g_{s} = \dfrac{D_{H_{2}O}}{D_{CO_{2}}} \dfrac{A_{n}}{c_{a}- c_{i}}
$$
```
where $\dfrac{D_{H_{2}O}}{D_{CO_{2}}} \approx 1.6$ is treated as a constant (`Drel`) in our model. 
## Citations
Chen, J.-L., Reynolds, J. F., Harley, P. C. & Tenhunen, J. D. Coordination theory of leaf nitrogen distribution in a canopy. Oecologia 93, 63–69 (1993)

Mengoli, G., Agustí-Panareda, A., Boussetta, S., Harrison, S. P., Trotta, C., & Prentice, I. C. (2022). Ecosystem photosynthesis in land-surface models: A first-principles approach incorporating acclimation. _Journal of Advances in Modeling Earth Systems_, 14, e2021MS002767. [https://doi.org/10.1029/2021MS002767](https://doi.org/10.1029/2021MS002767)

Stocker, B. D., Wang, H., Smith, N. G., Harrison, S. P., Keenan, T. F., Sandoval, D., Davis, T., & Prentice, I. C. (2020). P-model v1.0: An optimality-based light use efficiency model for simulating ecosystem gross primary production. _Geoscientific Model Development_, _13_(3), 1545–1581. [https://doi.org/10.5194/gmd-13-1545-2020](https://doi.org/10.5194/gmd-13-1545-2020)

Wang, H., Prentice, I. C., Keenan, T. F., Davis, T. W., Wright, I. J., Cornwell, W. K., Evans, B. J., & Peng, C. (2017). Towards a universal model for carbon dioxide uptake by plants. Nature Plants, 3(9), 734–741. https://doi.org/10.1038/s41477-017-0006-8
