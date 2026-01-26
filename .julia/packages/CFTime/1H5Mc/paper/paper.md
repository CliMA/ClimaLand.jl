---
title: 'CFTime.jl: a Julia package for representing time following the Climate and Forecast conventions'
tags:
  - julia
  - climate-and-forecast-conventions
  - oceanography
  - meteorology
  - earth-observation
  - climatology
  - netcdf
authors:
  - name: Alexander Barth
    orcid: 0000-0003-2952-5997
    affiliation: 1
affiliations:
 - name: GHER, University of Liège, Liège, Belgium
   index: 1
date: 13 January 2024
bibliography: paper.bib
---

# Summary


The Climate and Forecast (CF) conventions are a metadata standard for Earth data [@Eaton2024] and are mainly used in oceanography and meteorology.
The CF conventions were originally proposed for the NetCDF storage format, but they are also increasingly used with other formats like Zarr [@OGC_Zarr] and GRIB ([`GRIBDatasets.jl`](https://github.com/JuliaGeo/GRIBDatasets.jl)).

Since the initial release of the CF conventions [@Eaton2003], the encoding of time instances has been standardized. The Julia package `CFTime` implements various calendars that have been standardized in the frame of these conventions. It also supports some arithmetic operations for example, computing the duration between two time instances, ordering time instances and creating time ranges. The time origin and resolution are flexible ranging from days to attoseconds.


# Statement of need

In many Earth science disciplines and beyond, expressing a time instance and a duration is essential. The CF conventions provide a rich and flexible
framework for handling time, equally applicable to observations and model data. To our knowledge, `CFTime.jl` is the only package in the Julia ecosystem that implements the time structures standardized by the CF conventions. While almost all datasets used in Earth science use dates after the year 1582, some datasets or software systems  [e.g. @Octave; @SeaDataNet_format] use a time origin before this date, which makes it necessary to handle the transition from the Julian to the Gregorian calendar.
Some users also expressed the need for microseconds and nanoseconds as time resolutions even if they are rarely used in typical Earth science applications (`CFTime` github [issue 18](https://github.com/JuliaGeo/CFTime.jl/issues/18) and `NCDatasets` github [issue 248](https://github.com/JuliaGeo/NCDatasets.jl/issues/248)).

As of 18 September 2025, 132 Julia packages depend directly or indirectly on `CFTime` (excluding optional dependencies). For example, `CFTime` is used by numerical models, such as `ClimaOcean.jl`, a framework for realistic ocean and coupled sea-ice simulations based on `Oceananigans.jl` [@OceananigansJOSS], the hydrological modeling package `Wflow.jl` [@vanVerseveld2024] and `AIBECS.jl`, a modeling framework for global marine biogeochemical cycles [@Pasquier2022].

Several data-related packages also make direct or indirect use of `CFTime`, such as the NetCDF manipulation package `NCDatasets.jl` [@Barth2024], the gridded data processing package `YAXArrays.jl` [@Gans2023] and packages for handling in-situ data from various observing platforms (`OceanRobots.jl` [@Forget2024] and `ArgoData.jl` [@Forget2025]).

# State of the field

The API of `CFTime` was highly influenced by Julia's `Dates` module from the standard library. The `Dates` module implements the `DateTime` structure representing a time instance in the proleptic Gregorian calendar following the ISO 8601 standard. The time instances are encoded using a 64-bit integer with millisecond precision and 31 December 1 BC 00:00:00 as the starting date. In the Python ecosystem, the calendars from the CF conventions are implemented by the `cftime` package. The Python package implements the time to microsecond accuracy (as of version 1.6.4). A timestamp is represented internally by storing separately the year, month, day, hour, minute, second and microsecond. This reduces the risks of integer overflows at the expense of memory requirements and complexity of performing arithmetic operations.

# Installation

`CFTime` supports Julia 1.6, and later and can be installed with the Julia package manager using the following command:

```julia
using Pkg
Pkg.add("CFTime")
```
`CFTime` is a pure Julia package and currently depends only on the modules `Dates` and `Printf`, which are part of Julia’s standard library.

# Features

In the context of the CF conventions, a time instance is represented as a time offset measured from a time origin (in UTC): for example, the value 86400 with units "seconds since 1970-01-01 00:00:00" is 2 January 1970, 00:00:00 UTC. The units of the time offset and the time origin are stored in the `units` attribute of the time variable.

The `calendar` attribute of a NetCDF or Zarr time variable defines how the time offset and units are interpreted to derive the calendar year, month, day, hour, and so on.
The CF conventions define several calendar types, including:

| Calendar                | Type                         | Explanation |
| ----------------------- | ---------------------------- | ---------------------------- |
| `standard`, `gregorian` | `DateTimeStandard`           | the Gregorian calendar after 15 October 1582 and the Julian calendar before  |
| `proleptic_gregorian`   | `DateTimeProlepticGregorian` | the Gregorian calendar applied to all dates |
| `julian`                | `DateTimeJulian`             | the Julian calendar applied to all dates |
| `noleap`, `365_day`     | `DateTimeNoLeap`             | calendar without leap years |
| `all_leap`, `366_day`   | `DateTimeAllLeap`            | calendar with only leap years |
| `360_day`               | `DateTime360Day`             | calendar assuming that all months have 30 days |

`CFTime` is based on the Meeus' algorithm [@Meeus98] for the Gregorian and Julian calendars, with two adaptations:

* The original algorithm is based on floating-point arithmetic. The algorithm in `CFTime` is implemented using integer arithmetic, which is more efficient.
Additionally, underflows and overflows are easier to predict and handle with integer arithmetic.
* The Meeus' algorithm has been extended to dates prior to 100 BC.

The following is a list of the main features of `CFTime`:

* Basic arithmetic, such as subtracting two time instances to compute their duration, or adding a duration to a time instance.
* Support for a wide range of time resolutions, from days down to attoseconds, for feature parity with NumPy's `datetime64` type [@harris2020array; @numpy].
* Support for arbitrary time origins. Since the time origin for NumPy's `datetime64` type is fixed to 1 January 1970 at 00:00, the usefulness of some time units is limited. As an extreme example, with attoseconds, only a time span of ±9.2 s around the time origin can be represented since a 64-bit integer is used internally.
* By default, the time counter is a 64-bit integer, but other integer types or floating-point types can be used.

# Acknowledgements

I thank [all contributors](https://github.com/JuliaGeo/CFTime.jl/graphs/contributors) to this package, in particular Martijn Visser, Fabian Gans, Haakon Ludvig Langeland Ervik, Rafael Schouten and Yeesian Ng. I also acknowledge Jeff Whitaker for Python's [cftime](https://github.com/Unidata/cftime) which has helped the development of this package by providing reference values and a reference implementation for tests.

# Funding
The author acknowledges the F.R.S.-FNRS (Fonds de la Recherche Scientifique de Belgique) for funding his position. This work was partly performed with funding from the Blue-Cloud 2026 project under the Horizon Europe programme, Grant Agreement No. 101094227.

# References
