# Canopy Turbulent Fluxes

This page documents how turbulent fluxes are computed in the ClimaLand Canopy Model.

---

## Leaf properties

The leaf temperature is modelled in the Canopy Energy component. The leaf specific humidity is determined as the
saturated value at that temperature. The stomatal conductance at the canopy level is computed by the Canopy
Stomatal Conductance component. The problem is then to parameterize the air temperature and specific humidity
at the top of the canopy.
---

## Canopy top properties
We assume a flux balance between the sensible and vapor fluxes between the atmosphere and the canopy top, and
between the canopy top and the canopy leaves. For example,

```math
H  = -\rho_{\rm{atmos}} \frac{T_{\rm{atmos}}- T_{\rm{canopy, air}}}{r_{a,e}} = -\rho_{\rm{atmos}} \frac{T_{\rm{canopy, air}} - T_{\rm{leaves}}{r_{\rm{canopy}}},
```
and

```math
T  = -\rho_{\rm{atmos}} \frac{q_{\rm{atmos}}- q_{\rm{canopy, air}}}{r_{a,e}} = -\rho_{\rm{atmos}} \frac{q_{\rm{canopy, air}} - q_{\rm{leaves}}{r_{\rm{canopy}}+ r_{\rm{stomata}}}.
```

We can then solve for the canopy air properties:
```math
T_{\rm{canopy, air}} = ,
```
and

```math
q_{\rm{canopy, air}} = .
```
---

## Canopy Fluxes

Plugging the canopy air properties into the flux equations let's us write the fluxes between the canopy and the atmosphere as:
```math
H  = -\rho_{\rm{atmos}} \frac{T_{\rm{atmos}} - T_{\rm{leaves}}{r_{a,e}+ r_{\rm{canopy}}},
```
and

```math
T  = -\rho_{\rm{atmos}} \frac{q_{\rm{atmos}} - q_{\rm{leaves}}{r_{a,e} + r_{\rm{canopy}}+ r_{\rm{stomata}}}.
```

## Roughness lengths and displacement height


## Specifying canopy resistance

## Changing parameterizations
---
