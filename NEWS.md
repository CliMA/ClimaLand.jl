ClimaLand.jl Release Notes
========================

main
--------

- Moved neural snow to extension, reducing latency by a factor 2. PR
  [#567](https://github.com/CliMA/ClimaLand.jl/pull/567)

v0.11.1
-------
- ![][badge-✨feature] Add option to profile albedo job. PR
  [#505](https://github.com/CliMA/ClimaLand.jl/pull/551)

v0.11.0
-------
- Update SurfaceFluxes compat to include 0.11: PR [#548](https://github.com/CliMA/ClimaLand.jl/pull/548)
- Started changelog: PR [#547](https://github.com/CliMA/ClimaLand.jl/pull/547)
- ![][badge-🐛bugfix] Use ClimaCore v0.13.2, which includes remapping allocation fix for CPU case: PR [#546](https://github.com/CliMA/ClimaLand.jl/pull/546)
- ![][badge-💥breaking] Refactor BucketModelParameters, move some parameters to ClimaParams: PR [#507](https://github.com/CliMA/ClimaLand.jl/pull/507)
- ![][badge-✨feature] Add TOPMODEL runoff parameterization: PR [#511](https://github.com/CliMA/ClimaLand.jl/pull/511), Issue [#266](https://github.com/CliMA/ClimaLand.jl/issues/266)
- Add web dashboard for Fluxnet simulations: PR [#478](https://github.com/CliMA/ClimaLand.jl/pull/478)
- Use depot for buildkite pipeline: PR [#537](https://github.com/CliMA/ClimaLand.jl/pull/537)
- ![][badge-✨feature] Refactor EnergyHydrologyParameters, move some parameters to ClimaParams: PR [#505](https://github.com/CliMA/ClimaLand.jl/pull/505)
- Update to ClimaCore v0.13 (note: introduced CPU global bucket bug): PR [#536](https://github.com/CliMA/ClimaLand.jl/pull/536)
- Add ClimaLand logos: PR [#525](https://github.com/CliMA/ClimaLand.jl/pull/525)
- ![][badge-✨feature] Add infrastructure to use PFTs, add experiment running Ozark with PFTs: PR [#493](https://github.com/CliMA/ClimaLand.jl/pull/493)
- ![][badge-💥breaking] Refactor boundary condition types to unify types/structure: PR [#535](https://github.com/CliMA/ClimaLand.jl/pull/535)
- ![][badge-💥breaking] Update to ClimaParams v0.10, free Insolation and Thermodynamics, update SurfaceFluxes: PR [#508](https://github.com/CliMA/ClimaLand.jl/pull/508)
- ![][badge-✨feature] Add PrescribedPrecipitation driver to be used in global simulations: PR [#533](https://github.com/CliMA/ClimaLand.jl/pull/533)
- Switch to use new-central, pin Insolation (0.9.1) and Thermodynamics (0.12.3) versions: PR [#526](https://github.com/CliMA/ClimaLand.jl/pull/526)
- ![][badge-💥breaking] Refactor albedo parameterizations to have two albedo models, PrescribedBaregroundAlbedo and PrescribedSurfaceAlbedo: PR [#513](https://github.com/CliMA/ClimaLand.jl/pull/513)
- Add documentation of standalone models: PR [#474](https://github.com/CliMA/ClimaLand.jl/pull/474)

<!--
Contributors are welcome to begin the description of changelog items with badge(s) below. Here is a brief description of when to use badges for a particular pull request / set of changes:
 - 💥breaking - breaking changes. For example: removing deprecated functions/types, removing support for functionality, API changes, breaking changes in compats.
 - ✨feature - new feature added. For example: adding support for a new parameterization.
 - 🐛bugfix - bugfix. For example: fixing incorrect logic, resulting in incorrect results, or fixing code that otherwise might give a `MethodError`.
 - 🔥behavioralΔ - behavioral changes. For example: a new model is used, yielding more accurate results.
 - 🤖precisionΔ - machine-precision changes. For example, swapping the order of summed arguments can result in machine-precision changes.
 - 🚀performance - performance improvements. For example: improving type inference, reducing allocations, or code hoisting.
-->

[badge-💥breaking]: https://img.shields.io/badge/💥BREAKING-red.svg
[badge-✨feature]: https://img.shields.io/badge/feature/enhancement-blue.svg
[badge-🐛bugfix]: https://img.shields.io/badge/🐛bugfix-purple.svg
[badge-🔥behavioralΔ]: https://img.shields.io/badge/🔥behavioralΔ-orange.svg
[badge-🤖precisionΔ]: https://img.shields.io/badge/🤖precisionΔ-black.svg
[badge-🚀performance]: https://img.shields.io/badge/🚀performance-green.svg
