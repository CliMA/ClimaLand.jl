import ClimaUtilities.TimeVaryingInputs:
    TimeVaryingInput,
    AbstractTimeVaryingInput,
    LinearInterpolation,
    PeriodicCalendar
import Interpolations: Constant
import ClimaUtilities.Regridders: InterpolationsRegridder
import ClimaUtilities.FileReaders: NCFileReader, read
import ClimaUtilities.TimeManager: ITime, date
using Thermodynamics
using ClimaCore
using Dates
using Insolation
using SurfaceFluxes
import SurfaceFluxes.Parameters as SFP
using StaticArrays
import ..Parameters as LP
import ClimaLand

import ClimaComms
ClimaComms.@import_required_backends
context = ClimaComms.context()

start_date = DateTime("2008-03-01")
stop_date = DateTime("2008-06-01")

era5_ncdata_path = ClimaLand.Artifacts.find_era5_year_paths(
            start_date,
            stop_date;
            context,
        )
