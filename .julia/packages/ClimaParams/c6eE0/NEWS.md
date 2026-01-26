ClimaParams.jl Release Notes
========================
main
--------

v1.0.1
--------
- Add parameters for cloud top height diagnostic ([#251](https://github.com/CliMA/ClimaParams.jl/pull/251))

v1.0.0
--------
- Breaking: Removed the AbstractTOMLDict type ([#247](https://github.com/CliMA/ClimaParams.jl/pull/247))
- Added direct indexing `td["gravitational_acceleration"] ([#239](https://github.com/CliMA/ClimaParams.jl/pull/239))

v0.12.1
--------

- Log parameter when getindex is called ([#243](https://github.com/CliMA/ClimaParams.jl/pull/243))
- Change epoch_time to a DateTime ([#248](https://github.com/CliMA/ClimaParams.jl/pull/248))


v0.12.0
--------
- Change some parameter names ([#242](https://github.com/CliMA/ClimaParams.jl/pull/242))
- Change epoch_time to string ([#246](https://github.com/CliMA/ClimaParams.jl/pull/246))

v0.11.0
--------
- Update the specific heat of ice

v0.10.35
--------
- Add wet growth timescale parameter for P3 scheme ([#236](https://github.com/CliMA/ClimaParams.jl/pull/236))

v0.10.34
--------
- Add parameters for P3 scheme ([#235](https://github.com/CliMA/ClimaParams.jl/pull/235))
- Remove duplicate parameters for SB2006 ventilation parameterization ([#234](https://github.com/CliMA/ClimaParams.jl/pull/234))

v0.10.31
--------
- Add parameters for detrainment ramp and ustar^3 surface tke flux formulation. ([#231](https://github.com/CliMA/ClimaParams.jl/pull/231))

v0.10.30
--------
- Update the value of Prandtl max number ([#230](https://github.com/CliMA/ClimaParams.jl/pull/230))

v0.10.29
--------
- Add maximum value for Prandtl number ([#229](https://github.com/CliMA/ClimaParams.jl/pull/229))

v0.10.28
--------
- Add minimum value for inverse lambda parameter for 1-moment cloud microphysics ([#228](https://github.com/CliMA/ClimaParams.jl/pull/228))

v0.10.27
--------
- Add one more parameter for data driven number concnetration parameterization ([#227](https://github.com/CliMA/ClimaParams.jl/pull/227))

v0.10.26
--------
- Consolidate GCM-driven SCM relaxation parameters from separate shallow/deep versions into single parameters ([#225](https://github.com/CliMA/ClimaParams.jl/pull/225))

v0.10.25
--------
- Add parameter vectors for data-driven mixing length and perturbation pressure closures ([#220](https://github.com/CliMA/ClimaParams.jl/pull/220))
- Adds 8 timescale parameters for relaxing forced single column model profiles ([#210](https://github.com/CliMA/ClimaParams.jl/pull/210))

v0.10.24
--------
- Add parameters for data driven number concnetration parameterization ([#224](https://github.com/CliMA/ClimaParams.jl/pull/224))

v0.10.23
--------
- Add parameters for scaling factors for tracer hyperdiffusion and vertical diffusion. ([#223](https://github.com/CliMA/ClimaParams.jl/pull/223))

v0.10.22
--------
- Rename Middle Eastern Dust and add Asian Dust ABDINM. ([#221](https://github.com/CliMA/ClimaParams.jl/pull/221))

v0.10.21
--------
- Add parameter for constant slope parameterization in CloudMicrophysics.jl P3 scheme. ([#218](https://github.com/CliMA/ClimaParams.jl/pull/218))

v0.10.20
--------
- Add parameters for ABDINM and ABIFM for Asian, Saharan, Middle Eastern, and generic Dust ([#216](https://github.com/CliMA/ClimaParams.jl/pull/216))

v0.10.19
--------
- Add parameters for cloud effective radius diagnostics ([#215](https://github.com/CliMA/ClimaParams.jl/pull/215))

v0.10.18
--------
- Add parameters for ABDINM and ABIFM for Illite and Arizona Test Dust ([#214](https://github.com/CliMA/ClimaParams.jl/pull/214))

v0.10.17
--------
- Add parameters for DecayWithHeightDiffusion scheme ([#212](https://github.com/CliMA/ClimaParams.jl/pull/213))

v0.10.16
--------
- Add apparent snow and cloud ice density based on doi:10.1029/2020JD034157 ([#212](https://github.com/CliMA/ClimaParams.jl/pull/212))

v0.10.15
--------
- Add snow aspect ratio coefficients for Chen et al 2022 terminal velocity parameterization ([#211](https://github.com/CliMA/ClimaParams.jl/pull/211))

v0.10.14
------
- Remove all dependencies ([#208](https://github.com/CliMA/ClimaParams.jl/pull/208))

v0.10.13
------
- Remove mol co2 to kg C factor AutotrophicResp, add kg C to mol CO2 factor Heterotrophic Resp ([#205](https://github.com/CliMA/ClimaParams.jl/pull/205))

v0.10.12
------
- Add diagnostic covariance coeff, change turbulent entrainment parameter vec default  ([#204](https://github.com/CliMA/ClimaParams.jl/pull/204))

v0.10.11
------
- Add data-driven entrainment parameter vector ([#202](https://github.com/CliMA/ClimaParams.jl/pull/202))
- Add data-driven turbulent entrainment parameter vector ([#203](https://github.com/CliMA/ClimaParams.jl/pull/203))

v0.10.10
------
- Add parameters for surface wind gustiness ([#201](https://github.com/CliMA/ClimaParams.jl/pull/201))
- Update default values for edmf parameters ([#200](https://github.com/CliMA/ClimaParams.jl/pull/200))

v0.10.9
------
- Add parameters for RCEMIP surface temperature distribution ([#199](https://github.com/CliMA/ClimaParams.jl/pull/199))

v0.10.8
------
- Add parameters for minimum and maximum temperature for the radiation lookup table ([#198](https://github.com/CliMA/ClimaParams.jl/pull/198))
- Add example of Oceananigans single column model using ClimaParams ([#151](https://github.com/CliMA/ClimaParams.jl/pull/151))

v0.10.7
------
- Add CloudMicrophysics P3 parameters for sink terms ([#192](https://github.com/CliMA/ClimaParams.jl/pull/192))

v0.10.6
------
- Add 2-moment parameters for CloudMicrophysics ([#191](https://github.com/CliMA/ClimaParams.jl/pull/191))

v0.10.5
------
- Add linear fit parameters for homogeneous ice nucleation ([#366](https://github.com/CliMA/CloudMicrophysics.jl/pull/366))

v0.10.4
------
- Add albedo parameters ([#188](https://github.com/CliMA/ClimaParams.jl/pull/188))

v0.10.3
------
- Add Chen Terminal Velocity Table B5 coefficients ([#186](https://github.com/CliMA/ClimaParams.jl/pull/186))
- Add a detrainment parameter ([#187](https://github.com/CliMA/ClimaParams.jl/pull/187))

v0.10.2
-------
- Add parameters for Snow, MedlynConductance, and BeerLambert ([#164](https://github.com/CliMA/ClimaParams.jl/pull/164))
- Add SoilCO2Model Land parameters ([#179](https://github.com/CliMA/ClimaParams.jl/pull/179))

v0.10.1
-------
- Added Bucket model parameters ([#183](https://github.com/CliMA/ClimaParams.jl/pull/183))
- Added Energy Hydrology parameters ([#180](https://github.com/CliMA/ClimaParams.jl/pull/180))

v0.10.0
-------
- Renamed to ClimaParams ([#184](https://github.com/CliMA/ClimaParams.jl/pull/184))

v0.9.0
-------
- Started changelog
- Allow NamedTuples to be used as name maps ([#158](https://github.com/CliMA/ClimaParams.jl/pull/158))
- Update default value for `alpha_rayleigh_uh` ([#160](https://github.com/CliMA/ClimaParams.jl/pull/160))
- Add parameters for water based deposition nucleation for kaolinite, feldspar, and ferrihydrate ([#161](https://github.com/CliMA/ClimaParams.jl/pull/161))
- Fix typos in deposition nucleation parameters ([#162](https://github.com/CliMA/ClimaParams.jl/pull/162))
- Add parameters for the P3 scheme ([#163](https://github.com/CliMA/ClimaParams.jl/pull/163))
- Add autotrophic respiration parameters ([#165](https://github.com/CliMA/ClimaParams.jl/pull/165))
- Remove default type for TOML parsing ([#166](https://github.com/CliMA/ClimaParams.jl/pull/166))
- Replace and add additional ARG2000 parameters ([#130](https://github.com/CliMA/ClimaParams.jl/pull/130))
- Add `T_init_min` for thermodynamics saturation adjustment, changes T_min to 1 Kelvin 🧊 ([#171](https://github.com/CliMA/ClimaParams.jl/pull/171))
- Fix typos and group some parameters together ([#168](https://github.com/CliMA/ClimaParams.jl/pull/168))
- Add Frostenberg et al (2023) parameters ([#174](https://github.com/CliMA/ClimaParams.jl/pull/174))
