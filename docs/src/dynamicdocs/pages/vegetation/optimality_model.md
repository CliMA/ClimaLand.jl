# Optimality Model

Photosynthetic coordination theory, originally proposed by Von Caemmerer & Farquhar (1981), provides an approach to predict dynamic responses of photosynthetic capacity to environmental constraints. It primarily focuses on how leaf nitrogen (N) affects the photosynthetic capacity. Photosynthetic capacity varies both among plant types and over time and space, and a major determinant of photosynthetic capacity is the maximum rate of Rubisco carboxylation ($V_\text{cmax}$). 

In this optimality model, Smith et al. (2019) assumes that plants are able to acquire the N necessary to build leaves that can photosynthesize at the fastest possible rate given light availability and biophysical constraints. The Vcmax model estimates $V_\text{cmax}$ and $J_\text{max}$ as a function of environmental variables as follows:

```math
\begin{equation}
    V_\text{cmax}^* = \varphi I \left(\frac{m}{m_c}\right)\left(\frac{\overline{\omega}^*}{8\theta}\right),
\end{equation}
```
where
```math
\begin{equation}
    \overline{\omega}^* = 1 + \overline{\omega} - \sqrt{(1 + \overline{\omega})^2 - 4\theta\overline{\omega}},
\end{equation}
```
and
```math
\begin{equation}
    \overline{\omega} = -(1 - 2\theta) + \sqrt{(1 - \theta)\left(\frac{1}{\frac{4c}{m}\left(1 - \theta\frac{4c}{m}\right)} - 4\theta\right)},
\end{equation}
```
and
```math
\begin{equation}
    c = \frac{m}{8\theta}\left(1 - \frac{\varphi I + J_\text{max} - 2\theta\varphi I}{\sqrt{(\varphi I + J_\text{max})^2 - 4\theta\varphi I J_\text{max}}}\right)
\end{equation}
```
and
```math
\begin{equation}
    J_\text{max} = \varphi I \overline{\omega}
\end{equation}
```
```math
\begin{equation}
    m = \frac{C'_i - \Gamma^*}{C'_i + 2\Gamma^*}
\end{equation}
```
```math
\begin{equation}
    C'_i = \Gamma^* + (C_a - \Gamma^*)\frac{\xi}{\xi + \sqrt{D_g}}
\end{equation}
```
```math
\begin{equation}
    \xi = \sqrt{\beta \frac{K + \Gamma^*}{1.6\eta^*}}
\end{equation}
```
```math
\begin{equation}
    K = K_c\left(1 + \frac{O_i}{K_o}\right)
\end{equation}
```
```math
\begin{equation}
    m_c = \frac{C'_i - \Gamma^*}{C'_i + K}
\end{equation}
```

$\Gamma^*$ is the CO$_2$ compensation point in the absence of mitochondrial respiration
```math
\begin{equation}
    \Gamma^* = \Gamma^*_0 f(T, \Delta H_a) p/p_0
\end{equation}
```
where $\Gamma^*_0 = 4.332$ Pa, $p$ is the atmospheric pressure, $p_0 = 101325$ Pa, and $\Delta H_a = 37830$ J/mol.

The model has the following parameters:
- ``\varphi`` is the realized quantum yield of photosynthetic electron transport (dimensionless). Estimated at 0.257.
- ``\theta`` is the curvature of the light response curve (dimensionless). Estimated at 0.85.
- ``\beta`` is the ratio of the carbon cost of maintaining photosynthetic proteins to the carbon cost of maintaining a transpiration stream (dimensionless). Estimated at 146.

For Smith et al. (2019) Vcmax model:
- Altitude
- ``D_g`` is the vapor pressure deficit (VPD) at altitude
- ``C_a`` is the CO$_2$ partial pressure
- ``I`` is the incident photosynthetically active photon flux (PAR)
- ``T`` is the temperature
