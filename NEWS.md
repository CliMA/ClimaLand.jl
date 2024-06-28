ClimaLand.jl Release Notes
========================

main
--------
<<<<<<< HEAD
Improve Performance of soil canopy model
- PR [#666]
- PR [#677] Remove norm condition
Add base global soil canopy run to benchmark and experiments
- PR [#591] adds benchmark run
- PR [#669] adds topmodel runoff
Soil implicit solver
- PR [#678] updates Jacobian for soil energy
Misc
- PR [#685] Renames jacobian update functions to align with our convention
            for tendency functions
=======

v0.12.5
--------
- Turn of dss
  PR[#650](https://github.com/CliMA/ClimaLand.jl/pull/650)
- Add ClimaCoupler downstream test
  PR [#680](https://github.com/CliMA/ClimaLand.jl/pull/680)
- Adds ClimaLand.Diagnostics
  PR [#628](https://github.com/CliMA/ClimaLand.jl/pull/628)
- Improve Performance of soil canopy model
  PR [#666](https://github.com/CliMA/ClimaLand.jl/pull/666)
  PR [#677](https://github.com/CliMA/ClimaLand.jl/pull/677)
- Add base global soil canopy run to benchmark and experiments
  PR [#591](https://github.com/CliMA/ClimaLand.jl/pull/591)
  PR [#669](https://github.com/CliMA/ClimaLand.jl/pull/669)

>>>>>>> 668f23f0 (tag v0.12.5)
v0.12.4
--------
- Fix various canopy flux bugs
  PR [#641](https://github.com/CliMA/ClimaLand.jl/pull/641)
- Adds Richards to benchmark
  PR [#648](https://github.com/CliMA/ClimaLand.jl/pull/648)
- Reduce allocations in update runoffs 
  PR [#664](https://github.com/CliMA/ClimaLand.jl/pull/664)
- Store z in cache instead of model
  PR [#658](https://github.com/CliMA/ClimaLand.jl/pull/658)

v0.12.3
--------
- Add benchmark pipeline 
  PR [#592](https://github.com/CliMA/ClimaLand.jl/pull/592)
  PR [#642](https://github.com/CliMA/ClimaLand.jl/pull/642)
- Removed ArtifactWrappers
  PR [#627](https://github.com/CliMA/ClimaLand.jl/pull/627)
  PR [#540](https://github.com/CliMA/ClimaLand.jl/pull/640)
- Run RichardsModel on GPU
  PR [#638](https://github.com/CliMA/ClimaLand.jl/pull/638)

v0.12.2
--------
- Update implicit solver interface 
  PR [#542](https://github.com/CliMA/ClimaLand.jl/pull/542)
- Add ClimaLand.Artifacts
  PR [#624](https://github.com/CliMA/ClimaLand.jl/pull/624)

v0.12.1
--------
- Add `regridder_type` option to albedo constructors.
  PR [#662](https://github.com/CliMA/ClimaLand.jl/pull/622)
- ![][badge-✨feature] Add simple snow model.
  PR [#147](https://github.com/CliMA/ClimaLand.jl/pull/147)
- ![][badge-✨feature] Implemented a simple sublimation model, which is now
  included in the soil model. PR [#373](https://github.com/CliMA/ClimaLand.jl/pull/373)

v0.12.0
--------
- Updated to ClimaComms 0.6 and ClimaCore 0.14. Now, `CUDA` and `MPI` are no
  longer automatically installed.

v0.11.2
--------
-  Add profiling of soil/canopy model to buildkite pipeline. This was the
   first global run of this model, and we also fixed dss! for tuple-valued fields,
   changed how we compute the zenith angle to allow for 2d fields, and changed
   instances of `sum` to a column integral.
   PR [#561](https://github.com/CliMA/ClimaLand.jl/pull/561)
- ![][badge-✨feature] Use
  [ClimaUtilities](https://github.com/CliMA/ClimaUtilities.jl) for `Space` and
  `Time` `VaryingInputs`. This adds support to non-conservative MPI/GPU
  compatible regridding of NetCDF input. PR
  [#560](https://github.com/CliMA/ClimaLand.jl/pull/560)
- Moved neural snow to extension, reducing latency by a factor 2. PR
  [#567](https://github.com/CliMA/ClimaLand.jl/pull/567)

v0.11.1
-------
- ![][badge-✨feature] Add option to profile albedo job. PR
  [#505](https://github.com/CliMA/ClimaLand.jl/pull/551)
- ![][badge-✨feature] ClimaLandSimulations: better plots (legend and style tweaks, added cumulative ET and P), unit tests, start and end time as an optional argument. PR [#494](https://github.com/CliMA/ClimaLand.jl/pull/494)
- ![][badge-🐛bugfix] ClimaLandSimulations: installation readme. PR [#494](https://github.com/CliMA/ClimaLand.jl/pull/494)

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
