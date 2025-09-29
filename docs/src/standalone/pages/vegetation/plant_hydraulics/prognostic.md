# Plant Hydraulics

Water loss during day-time transpiration drives plants to draw water from the soil by roots and transport it through the stem to leaves. The plant hydraulics code solves for the volumetric water content along the water flow path, in the stem and leaf ($\theta_{stem}$ and $\theta_{leaf}$). It allows for an arbitrary number of compartments, but for now we will start with a single stem and leaf compartment. 

The volume flux of water $q$ (m/s) between compartments is given by Darcy's law as
```math
\begin{align}
    q = - K(\psi) \frac{dh}{dz}
\end{align}
```
where $h = \psi+z$ is the head (in meters), and $K$ is the conductivity (units of m/s). We approximate the flux between two compartments, indexed with 2 at
height $z_2$ and indexed with 1 at height $z_1<z_2$, using finite difference as
```math
\begin{equation}
q = \approx -\frac{k_1(\psi_1) k_2(\psi_2)}{k_1(\psi_1)+k_2(\psi_2)} * \bigg[\frac{\psi_2 - \psi_1}{z_2 - z_1} +1\bigg].
\end{equation}
```
The change of water volume (m$^3),  V$, in the compartments is then
```math
\begin{align}
    \frac{d V_{w, stem}}{dt} = q_{roots}\sigma_{roots} - q_{stem}\sigma_{stem} \nonumber \\
    \frac{d V_{w, leaf}}{dt} = q_{stem}\sigma_{stem} - \tau \sigma_{leaf},
\end{align}
```
where $\tau$ is a transpiration volume flux per unit emitting area, and $\sigma$ is the total emitting/conducting area.

This currently holds for a single pathway. To convert to fluxes from an entire surface, we can multiply by the number of individuals $N$. We can make use of the fact that $N\sigma/A$, where $A$ is the area of the ground those $N$ individuals are occupying, is the area index for that plant type. Following CLM, we incorporate a root, stem, and leaf area index (RAI, SAI, LAI) in order to model fluxes across an entire grid cell. 

Then we have:
```math
\begin{align}
    \frac{d v_{stem}}{dt} = q_{roots}RAI - q_{stem}SAI \nonumber \\
    \frac{d v_{leaf}}{dt} = q_{stem}SAI - \tau LAI,
\end{align}
```
where $v$ now represents the volume of water in that compartment (of a bulk plant) per unit ground area.

We also need to convert from the variable $v$ to $\psi$, in order to compute root extraction with the soil.  To do so, we can convert $v$ to the volumetric water content, and from $\theta$ to $\psi$ using the retention curve. To convert, let the volume of water per area of compartment be $V_{w}$, and $H$ the typical ``length" of the compartment. Then e.g.,
```math
\begin{equation}
    \theta_{stem}=\frac{V_{w,stem}}{A_{ground}} \times \frac{A_{ground}}{A_{stem}} \times \frac{1}{H_{stem}} = \frac{v_{stem}}{H_{stem} \times SAI }.
\end{equation}
```
Substituting in the volumetric water content, we have
```math
\begin{align}
    \frac{d \theta_{stem}}{dt} &= \frac{q_{roots}RAI - q_{stem}SAI}{H_{stem} SAI} \nonumber \\
    \frac{d \theta_{leaf}}{dt} &= \frac{q_{stem}SAI - \tau LAI}{H_{leaf} LAI},
\end{align}
```

We can also account for the distribution of roots as a function of depth. A quantity that is modeled in plant hydraulic models is the root fraction $P(z)$, satisfying $\int P(z) dz = 1$. Instead of having a single root at one discrete location, we can distribute the root system over different depths using $P(z)$. The total flux from roots between $z$ and $z+dz$ is given by
```math
\begin{equation}
    dq_{roots}(z) = -P(z) q(z) dz.
\end{equation}
```
Then the net root flux for the plant system would sum over this
```math
\begin{equation}
    q_{roots}  = -\int_{z_{min}}^{z_{sfc}} P(z) q(z) dz ,
\end{equation}
```
where $z_{min}$ is the minimum soil layer of the simulation. 

The sink term of the soil is in terms of a volumetric fraction change, i.e. we need a volume of water per volume of soil per second. We can obtain this with
```math
\begin{equation}
    S(z) = -(RAI) dq_{roots}(z)/dz = (RAI) P(z) q(z).
\end{equation}
```

The sign change occurs in the expression for $S(z)$ because a positive value of $q_{roots}$ indicates flow from the soil to the plant. This is a sink term for the soil.

In order to close the set of equations, the user will have to specify $K(\psi)$ and a function $\psi(\theta)$. In our current implementation, we use the same parameter values for all compartments, a linear retention curve $\psi(\theta)$, and a Weibull permeability curve $K(\psi)$:
```math
\begin{equation}
    K(\psi)  = \begin{cases} K_{\rm sat}\exp{(\psi/\psi_{63})^c} & \psi <0 \\
K_{\rm sat} & \psi > 0
\end{cases}
\end{equation}
```
and
```math
\begin{equation}
    \psi(\theta)  = \frac{1}{a}(\theta/\nu - 1).
\end{equation}
```
The model needs the following parameters:

| Drivers | Symbol | Unit | Approximate Value |
| :---         |     :---:      |    :---:      |     :---:   |
| Linear retention curve slope | a | 1/m | 0.002 |
| Saturated conductivity | $K_{\rm sat}$ | m/s  | 1e-7 |
| Potential at 1/e loss in conductivity | $\p_{63}$ | m | -400 |
| Weibull exponential parameter | $c$ | - | 4 |
| Porosity | $\nu$ | - | 1e-4 |
