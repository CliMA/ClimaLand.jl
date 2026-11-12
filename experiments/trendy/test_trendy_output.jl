# Tests of experiments/trendy/trendy_output.jl on synthetic ClimaDiagnostics files.
#   julia --project=.buildkite experiments/trendy/test_trendy_output.jl
using Test
import NCDatasets
import Dates

include(joinpath(@__DIR__, "trendy_output.jl"))

"""
Write a monthly ClimaDiagnostics-like file `<name>_1M_average.nc` in `dir`, with
values `f(lon, lat, [z,] month index)` on land (|lat| < 60) and NaN elsewhere.
"""
function write_fake_diagnostic(
    dir,
    name,
    lon,
    lat,
    f;
    z = nothing,
    nmonths = 14,
)
    start_date = Dates.DateTime(2003)
    starts = [start_date + Dates.Month(m - 1) for m in 1:nmonths]
    stops = [start_date + Dates.Month(m) for m in 1:nmonths]
    seconds(d) = Dates.value(Dates.Second(d - start_date))
    NCDatasets.NCDataset(joinpath(dir, "$(name)_1M_average.nc"), "c") do ds
        NCDatasets.defVar(ds, "lon", lon, ("lon",))
        NCDatasets.defVar(ds, "lat", lat, ("lat",))
        NCDatasets.defVar(
            ds,
            "time",
            Float64.(seconds.(starts)),
            ("time",);
            attrib = ["units" => "seconds since 2003-01-01 00:00:00"],
        )
        NCDatasets.defVar(
            ds,
            "time_bnds",
            Float64.(vcat(seconds.(starts)', seconds.(stops)')),
            ("nv", "time");
            attrib = ["units" => "seconds since 2003-01-01 00:00:00"],
        )
        if isnothing(z)
            data = [
                abs(φ) < 60 ? f(λ, φ, t) : NaN for
                λ in lon, φ in lat, t in 1:nmonths
            ]
            NCDatasets.defVar(ds, name, data, ("lon", "lat", "time"))
        else
            NCDatasets.defVar(ds, "z", z, ("z",))
            data = [
                abs(φ) < 60 ? f(λ, φ, zk, t) : NaN for
                λ in lon, φ in lat, zk in z, t in 1:nmonths
            ]
            NCDatasets.defVar(ds, name, data, ("lon", "lat", "z", "time"))
        end
    end
end

@testset "TRENDY output writer" begin
    diagnostics_dir = mktempdir()
    output_root = mktempdir()
    # ClimaDiagnostics-like 2° grid, longitudes 0 to 358
    src_lon = collect(0.0:2.0:358.0)
    src_lat = collect(-89.0:2.0:89.0)
    z = [-0.05, -0.15, -0.5, -1.5] # layer centres; faces 0, -0.1, -0.325, -1, -2
    gpp_mol = 1e-5
    write_fake_diagnostic(
        diagnostics_dir,
        "gpp",
        src_lon,
        src_lat,
        (λ, φ, t) -> gpp_mol,
    )
    write_fake_diagnostic(
        diagnostics_dir,
        "ra",
        src_lon,
        src_lat,
        (λ, φ, t) -> 0.4gpp_mol,
    )
    write_fake_diagnostic(
        diagnostics_dir,
        "hr",
        src_lon,
        src_lat,
        (λ, φ, t) -> 0.5gpp_mol,
    )
    # nee = ra + hr - gpp: a sink of 0.1 gpp
    write_fake_diagnostic(
        diagnostics_dir,
        "nee",
        src_lon,
        src_lat,
        (λ, φ, t) -> -0.1gpp_mol,
    )
    write_fake_diagnostic(
        diagnostics_dir,
        "precip",
        src_lon,
        src_lat,
        (λ, φ, t) -> -3e-5,
    )
    write_fake_diagnostic(
        diagnostics_dir,
        "tr",
        src_lon,
        src_lat,
        (λ, φ, t) -> 1e-8,
    )
    write_fake_diagnostic(
        diagnostics_dir,
        "tair",
        src_lon,
        src_lat,
        (λ, φ, t) -> 280.0 + φ / 10,
    )
    write_fake_diagnostic(
        diagnostics_dir,
        "soc",
        src_lon,
        src_lat,
        (λ, φ, zk, t) -> 10.0;
        z,
    )
    write_fake_diagnostic(
        diagnostics_dir,
        "swc",
        src_lon,
        src_lat,
        (λ, φ, zk, t) -> 0.3;
        z,
    )
    write_fake_diagnostic(
        diagnostics_dir,
        "si",
        src_lon,
        src_lat,
        (λ, φ, zk, t) -> 0.1;
        z,
    )
    write_fake_diagnostic(
        diagnostics_dir,
        "tsoil",
        src_lon,
        src_lat,
        (λ, φ, zk, t) -> 270.0 - zk;
        z,
    )

    @test layer_thickness(z) ≈ [0.1, 0.225, 0.675, 1.0]
    @test layer_thickness(z; soil_depth = 3.0) ≈ [0.1, 0.225, 0.675, 2.0]
    @test layer_faces(z; soil_depth = 3.0) ≈ [0.0, -0.1, -0.325, -1.0, -3.0]
    @test top_layer_weights(z, 0.1) ≈ [1.0, 0.0, 0.0, 0.0]
    @test top_layer_weights(z, 0.2) ≈ [1.0, 0.1 / 0.225, 0.0, 0.0]
    # Land fraction: a 1° land cell centred on integer degrees covers a quarter of
    # each of the four 1° TRENDY cells it overlaps; totals are conserved
    tlon, tlat = trendy_grid(1.0)
    ilon = collect(-180.0:179.0)
    ilat = collect(-90.0:89.0)
    src = zeros(360, 180)
    src[181, 101] = 1.0 # centred at 0°E, 10°N
    frac = land_fraction_on_grid(src, ilon, ilat, tlon, tlat)
    @test sum(frac .> 0) == 4
    @test sum(frac .* grid_cell_area(tlon, tlat)) ≈
          EARTH_RADIUS^2 * deg2rad(1) * (sind(10.5) - sind(9.5)) rtol = 1e-10
    # The source grid has no row at 90°N, so the band 89.5–90°N is not covered
    full = land_fraction_on_grid(ones(360, 180), ilon, ilat, tlon, tlat)
    @test full[:, 1:179] ≈ ones(360, 179)
    @test full[:, 180] ≈ fill((sind(89.5) - sind(89)) / (1 - sind(89)), 360)
    lon, lat = trendy_grid(0.5)
    @test length(lon) == 720 && length(lat) == 360
    @test lon[1] == -179.75 && lat[end] == 89.75
    @test sum(grid_cell_area(trendy_grid(1.0)...)) ≈ 4π * EARTH_RADIUS^2 rtol =
        1e-10

    @test_throws ErrorException write_trendy_output(
        diagnostics_dir,
        output_root;
        model = "Clima_Land",
    )
    written = write_trendy_output(
        diagnostics_dir,
        output_root;
        model = "ClimaLand-OBS",
        simulation = "S2",
        resolution = 1.0,
        commit = "abc1234",
    )
    sim_dir = joinpath(output_root, "ClimaLand-OBS", "S2")
    for name in (
        "gpp",
        "ra",
        "npp",
        "rh",
        "nbp",
        "fLuc",
        "pr",
        "mrro",
        "tas",
        "cSoil",
        "mrso",
        "mrsos",
        "tsl",
        "msl",
        "oceanCoverFrac",
    )
        @test isfile(joinpath(sim_dir, "ClimaLand-OBS_S2_$(name).nc"))
    end
    @test !isfile(joinpath(sim_dir, "ClimaLand-OBS_S2_lai.nc"))
    @test isfile(
        joinpath(output_root, "ClimaLand-OBS", "README_ClimaLand-OBS.txt"),
    )

    NCDatasets.NCDataset(joinpath(sim_dir, "ClimaLand-OBS_S2_nbp.nc")) do ds
        @test NCDatasets.dimnames(ds["nbp"]) ==
              ("longitude", "latitude", "time")
        @test ds["longitude"][1] == -179.5 && ds["latitude"][1] == -89.5
        @test ds["nbp"].attrib["units"] == "kg m-2 s-1"
        @test ds["nbp"].attrib["_FillValue"] == -99999.0f0
        @test ds["nbp"].attrib["missing_value"] == -99999.0f0
        @test ds["time"].attrib["calendar"] == "365_day"
        @test ds.attrib["commit"] == "abc1234"
        raw = ds["nbp"].var[:, :, :]
        @test raw[1, 1, 1] == -99999.0f0 # ocean (lat < -60) is fill
        @test raw[100, 90, 1] ≈ 0.1gpp_mol * KG_C_PER_MOL_CO2 rtol = 1e-6
        # Mid-January and mid-February 2003 in the 365-day calendar
        @test ds["time"].var[1:2] == [365 * 303 + 15.5, 365 * 303 + 31 + 14]
        @test ds["time_bnds"].var[:, 1] == [365 * 303, 365 * 303 + 31]
        @test length(ds["time"]) == 14
    end
    NCDatasets.NCDataset(joinpath(sim_dir, "ClimaLand-OBS_S2_npp.nc")) do ds
        @test ds["npp"].var[100, 90, 1] ≈ 0.6gpp_mol * KG_C_PER_MOL_CO2 rtol =
            1e-6
    end
    NCDatasets.NCDataset(joinpath(sim_dir, "ClimaLand-OBS_S2_pr.nc")) do ds
        @test ds["pr"].var[100, 90, 1] ≈ 3e-5 rtol = 1e-6
    end
    NCDatasets.NCDataset(joinpath(sim_dir, "ClimaLand-OBS_S2_fLuc.nc")) do ds
        @test ds["fLuc"].var[100, 90, 1] == 0
        @test ds["fLuc"].var[1, 1, 1] == -99999.0f0
    end
    NCDatasets.NCDataset(joinpath(sim_dir, "ClimaLand-OBS_S2_cSoil.nc")) do ds
        @test ds["cSoil"].var[100, 90, 1] ≈ 10.0 * 2.0 rtol = 1e-6
    end
    NCDatasets.NCDataset(joinpath(sim_dir, "ClimaLand-OBS_S2_mrso.nc")) do ds
        @test ds["mrso"].var[100, 90, 1] ≈ (300.0 + 91.67) * 2.0 rtol = 1e-6
        @test ds["mrso"].attrib["soil_column_depth"] ≈ 2.0
    end
    NCDatasets.NCDataset(joinpath(sim_dir, "ClimaLand-OBS_S2_mrsos.nc")) do ds
        @test ds["mrsos"].var[100, 90, 1] ≈ (300.0 + 91.67) * 0.1 rtol = 1e-6
        @test ds["mrsos"].attrib["soil_depth"] ≈ 0.1
    end
    NCDatasets.NCDataset(joinpath(sim_dir, "ClimaLand-OBS_S2_msl.nc")) do ds
        @test ds["smlayer_bnds"][:, :] ≈
              [0.0 0.1 0.325 1.0; 0.1 0.325 1.0 2.0]
        @test sum(ds["msl"].var[100, 90, :, 1]) ≈ (300.0 + 91.67) * 2.0 rtol =
            1e-5
    end
    NCDatasets.NCDataset(joinpath(sim_dir, "ClimaLand-OBS_S2_tsl.nc")) do ds
        @test NCDatasets.dimnames(ds["tsl"]) ==
              ("longitude", "latitude", "stlayer", "time")
        @test ds["stlayer"][:] ≈ -z
        @test ds["tsl"].var[100, 90, 4, 1] ≈ 271.5
    end
    NCDatasets.NCDataset(joinpath(sim_dir, "ClimaLand-OBS_S2_tas.nc")) do ds
        # Bilinear from the 2° grid: at 0.5°N, between 280.0 + 0.1 * (-1, 1)
        @test ds["tas"].var[181, 91, 1] ≈ 280.05 rtol = 1e-6
    end
    NCDatasets.NCDataset(
        joinpath(sim_dir, "ClimaLand-OBS_S2_oceanCoverFrac.nc"),
    ) do ds
        @test NCDatasets.dimnames(ds["oceanCoverFrac"]) ==
              ("longitude", "latitude")
        @test ds["oceanCoverFrac"][1, 1] == 1 &&
              ds["oceanCoverFrac"][100, 90] == 0
    end

    # Zonal NBP: only the complete year 2003; uniform sink over the land, whose 2°
    # cells (centres up to ±59°) end at ±60°; the land fraction conserves that area
    zonal = readlines(
        joinpath(output_root, "ClimaLand-OBS", "ClimaLand-OBS_zonalNBP_S2.txt"),
    )
    @test zonal[1] == "Year\tGlobal\tNorth\tTropics\tSouth"
    @test length(zonal) == 2
    values = parse.(Float64, split(zonal[2], '\t')[2:end])
    land_area(φ1, φ2) = 2π * EARTH_RADIUS^2 * (sind(φ2) - sind(φ1))
    flux = 0.1gpp_mol * KG_C_PER_MOL_CO2 * 365 * 86400 * 1e-12
    @test values[1] ≈ flux * land_area(-60, 60) rtol = 1e-3
    @test values[2] ≈ flux * land_area(30, 60) rtol = 1e-3
    @test values[3] ≈ flux * land_area(-30, 30) rtol = 1e-3
    @test values[1] ≈ sum(values[2:4]) rtol = 1e-4

    # 0.5° output and restricted variables
    written = write_trendy_output(
        diagnostics_dir,
        output_root;
        model = "ClimaLand",
        simulation = "S3",
        resolution = 0.5,
        variables = ["gpp"],
    )
    NCDatasets.NCDataset(
        joinpath(output_root, "ClimaLand", "S3", "ClimaLand_S3_gpp.nc"),
    ) do ds
        @test length(ds["longitude"]) == 720 && length(ds["latitude"]) == 360
    end
    @test !isfile(
        joinpath(output_root, "ClimaLand", "S3", "ClimaLand_S3_nbp.nc"),
    )
    @test occursin("| `nbp` |", trendy_mapping_table())
end
