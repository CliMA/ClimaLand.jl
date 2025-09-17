# Canopy Structure (Prescribed)

This page documents the **prescribed canopy structure/biomass** component used by ClimaLand’s vegetation module. This module does **not** prognose vegetation; it simply exposes structural and biomass properties to other sub‑models (radiative transfer, stomata/photosynthesis, energy, hydraulics) as **given inputs**.

It corresponds to `src/standalone/Vegetation/biomass.jl`.

---

## What it contains

- **Structural state**
  - Total leaf area index (**LAI**), stem area index (SAI), and root area index (RAI).
  - Canopy height, Rooting depth.

---
