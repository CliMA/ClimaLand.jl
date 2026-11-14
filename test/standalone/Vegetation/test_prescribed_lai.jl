using Test
using Dates
import ClimaComms
ClimaComms.@import_required_backends
using ClimaCore
import ClimaUtilities.TimeManager: ITime
import ClimaUtilities.TimeVaryingInputs: evaluate!
using ClimaLand
using ClimaLand.Canopy: MidMonthInput, stamp_date

FT = Float32

@testset "MODIS LAI monthly means at mid-month" begin
    start_date = DateTime(2008, 3, 1)
    stop_date = DateTime(2009, 3, 1)
    domain = ClimaLand.Domains.global_domain(FT; nelements = (10, 5))
    surface_space = domain.space.surface
    lai = ClimaLand.Canopy.prescribed_lai_modis(
        surface_space,
        start_date,
        stop_date,
    )
    raw = ClimaLand.Canopy.prescribed_lai_modis(
        surface_space,
        start_date,
        stop_date;
        midmonth = false,
    )
    @test lai isa MidMonthInput
    @test !(raw isa MidMonthInput)

    # The files stamp month k of a year at January 1 + 30(k - 1) days
    @test lai.stamps[7] == DateTime(2008, 6, 29)
    @test lai.middles[7] == DateTime(2008, 7, 16, 12)
    @test lai.middles[2] == DateTime(2008, 2, 15)
    @test lai.middles[12] == DateTime(2008, 12, 16, 12)
    @test lai.middles[13] == DateTime(2009, 1, 16, 12)

    # A mid-month value is the monthly mean; halfway between two middles is
    # halfway between the two stamps
    @test stamp_date(lai, DateTime(2008, 7, 16, 12)) == DateTime(2008, 6, 29)
    halfway = lai.middles[7] + (lai.middles[8] - lai.middles[7]) ÷ 2
    @test stamp_date(lai, halfway) ==
          lai.stamps[7] + (lai.stamps[8] - lai.stamps[7]) ÷ 2
    @test stamp_date(lai, DateTime(2007, 12, 1)) == lai.stamps[begin]
    @test stamp_date(lai, DateTime(2010, 6, 1)) == lai.stamps[end]

    # Stamps already at mid-month are kept
    @test stamp_date(
        MidMonthInput(lai.input, lai.middles, lai.middles),
        halfway,
    ) == halfway

    value = ClimaCore.Fields.zeros(surface_space)
    expected = ClimaCore.Fields.zeros(surface_space)
    t = Float64(Dates.value(Second(DateTime(2008, 7, 16, 12) - start_date)))
    evaluate!(value, lai, t)
    evaluate!(
        expected,
        raw,
        Float64(Dates.value(Second(DateTime(2008, 6, 29) - start_date))),
    )
    @test Array(parent(value)) ≈ Array(parent(expected))
    evaluate!(value, lai, ITime(round(Int, t), Second(1), start_date))
    @test Array(parent(value)) ≈ Array(parent(expected))
end
