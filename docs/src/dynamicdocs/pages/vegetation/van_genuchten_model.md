# Plant Hydraulics
The plant hydraulics code solves for the volumetric water content in the stem and leaf ($\theta_{stem}$ and $\theta_{leaf}$). It allows for an arbitrary number of stem/leaf compartments, but for now we will start with a single stem and leaf compartment. 

## Van Genuchten Model
The volume flux of water $q$ (m/s) between compartments with centers at two heights, $z_1$ and $z_2$, is given by Darcy's law as
```math
\begin{align}
    q = -\int_{z_1}^{z_2} k(\psi) dh
\end{align}
```
where $h = \psi+z$ is the head (in meters), and $k$ is the conductance (units of 1/s). As this is the conductance unit that CLM uses, there should be data bases with this information. We approximate this using finite difference as\footnote{Double check this - the units of $k$ in our code are $m/s$.}
```math
\begin{equation}
q = -\int_{h_1}^{h_2} k(\psi) dh \approx -\frac{k_1(\psi_1) + k_2(\psi_2)}{2} * [(\psi_2 - \psi_1) + (z_2 - z_1)].
\end{equation}
```
In order to close the set of equations, the user will have to specify $k(\psi)$ and a function $\psi(\theta)$. In our current implementation, we use a van Genuchten relationship with the same parameters for all compartments, but differing values of $K_{sat}.$

The change of water volume (m$^3),  V$, in the compartments is then
```math
\begin{align}
    \frac{d V_{w, stem}}{dt} = q_{roots}\sigma_{roots} - q_{stem}\sigma_{stem} \nonumber \\
    \frac{d V_{w, leaf}}{dt} = q_{stem}\sigma_{stem} - \tau \sigma_{leaf},
\end{align}
```
where $\tau$ is a transpiration volume flux per unit emitting area, and $\sigma$ is the total emitting/conducting area\footnote{Note that these are actually the areas at the faces between compartments. In the code, we take the average of the cross section of the compartments to estimate this.}. 

This currently holds for a single plant. To convert to fluxes from an entire surface, we can multiply by the number of individuals $N$. We can make use of the fact that $N\sigma/A$, where $A$ is the area of the ground those $N$ individuals are occupying, is the area index for that plant type. Following CLM, we incorporate a root, stem, and leaf area index (RAI, SAI, LAI) in order to model fluxes across an entire grid cell. 

Then we have:
```math
\begin{align}
    \frac{d v_{stem}}{dt} = q_{roots}RAI - q_{stem}SAI \nonumber \\
    \frac{d v_{leaf}}{dt} = q_{stem}SAI - \tau LAI,
\end{align}
```
where $v$ now represents the volume of water in that compartment (of a bulk plant) per unit ground area.

We also need to convert from the variable $v$ to $\psi$, in order to compute root extraction with the soil.  To do so, we can convert $v$ to the volumetric water content, and from $\theta$ to $\psi$ using a van Genuchten relationship. To convert, let the volume of water per area of compartment be $V_{w,stem}$, and $H$ the typical ``length" of the compartment. Then
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
    dq_{roots}(z) = -P(z) dz \int_{h_{soil}(z)}^{h_{stem}} k(\psi) dh,
\end{equation}
```
so that the net flux for the plant system would sum over this
```math
\begin{equation}
    q_{roots}  = -\int_{z_{min}}^{z_{sfc}} \frac{dq_{roots}(z)}{dz}dz ,
\end{equation}
```
where $z_{min}$ is the minimum soil layer of the simulation. 

The sink term of the soil is in terms of a volumetric fraction change, i.e. we need a volume of water per volume of soil per second. We can obtain this with
```math
\begin{equation}
    S(z) = -(RAI) dq_{roots}(z)/dz = (RAI) P(z) \int_{h_{soil}(z)}^{h_{stem}} k(\psi) dh.
\end{equation}
```

The sign change occurs in the expression for $S(z)$ because a positive value of $q_{roots}$ indicates flow from the soil to the plant. This is a sink term for the soil.
