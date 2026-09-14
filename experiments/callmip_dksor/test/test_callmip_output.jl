#=
Format-compliance test for the DK-Sor CalLMIP Phase 1a NetCDF writer.

Builds tiny synthetic diagnostics, runs `write_callmip_nc` (reused from
write_callmip_netcdf.jl via its PROGRAM_FILE include-guard), and asserts the
output matches the official CalLMIP template
(template_output_file/CLASSIC.v2_Phase1b_Scen1_DK-Sor_Cal_Prior.nc) EXACTLY where
it must: variable names, units, dimension order (time,lat,lon), time axis
("days since 1996-12-31 00:00", values 1..6574, standard calendar), and the
required Phase 1a global attributes.

Run:  julia --project=.buildkite experiments/callmip_dksor/test/test_callmip_output.jl
=#
import ClimaComms; ClimaComms.@import_required_backends
using Test, JLD2, NCDatasets, Dates

include(joinpath(@__DIR__, "..", "write_callmip_netcdf.jl"))  # PROGRAM_FILE guard ⇒ no main run

# Expected spec, verified against the CalLMIP template .nc
const EXPECT_UNITS = Dict(
    "nep" => "kgC m-2 s-1", "hfls" => "W m-2", "hfss" => "W m-2",
    "gpp" => "kgC m-2 s-1", "reco" => "kgC m-2 s-1", "tran" => "kg m-2 s-1",
    "evspsblsoi" => "kg m-2 s-1", "hfg" => "W m-2", "ts" => "K",
    "mrso" => "kg m-2", "mrsos" => "kg m-2", "lai" => "m2 m-2",
    "cSoil" => "kg m-2", "cLiveBioAbove" => "kg m-2",
)

# ── synthetic diagnostics (3 days, 5 soil layers) ─────────────────────────────
dir = mktempdir()
dates = [Date(2009, 6, d) for d in 1:3]
sd = Dict{String,Vector{Float64}}(k => Float64.(1:3) for k in
     ("nee", "lhf", "shf", "gpp", "er", "trans", "ct", "lai", "cveg"))
cd_ = Dict("swc" => fill(0.3, 5, 3), "tsoil" => fill(285.0, 5, 3), "soc" => fill(2.0, 5, 3))
z = Float64[-0.05, -0.2, -0.5, -1.0, -2.0]
diag = joinpath(dir, "callmip_diagnostics.jld2")
jldsave(diag; dates = dates, surface_data = sd, column_data = cd_, z_soil = z)

ncf = joinpath(dir, "test_out.nc")
write_callmip_nc(diag, nothing, ncf, "Prior")

@testset "CalLMIP Phase 1a NetCDF format" begin
    NCDataset(ncf, "r") do ds
        @testset "dimensions + time axis" begin
            @test ds.dim["time"] == 6574          # ds.dim[name] returns the length
            @test ds.dim["lat"] == 1 && ds.dim["lon"] == 1
            t = ds["time"].var[:]
            @test t[1] == 1.0 && t[end] == 6574.0
            @test ds["time"].attrib["units"] == "days since 1996-12-31 00:00"
            @test ds["time"].attrib["calendar"] == "standard"
        end
        @testset "data variables: presence, units, dim order" begin
            for (v, u) in EXPECT_UNITS
                @test haskey(ds, v)
                @test ds[v].attrib["units"] == u
                # NCDatasets (col-major) reports the reverse of on-disk order, so
                # ("lon","lat","time") here ⇒ on-disk/CF (time,lat,lon) as required.
                @test NCDatasets.dimnames(ds[v]) == ("lon", "lat", "time")
                @test ds[v].attrib["_FillValue"] == 1.0e38
            end
        end
        @testset "required global attributes (Phase 1a)" begin
            @test ds.attrib["Model"] == "ClimaLand.v1"
            @test ds.attrib["CalLMIP_Phase"] == "Phase1a"
            @test ds.attrib["Calibration_Scenario"] == "Scen1"
            @test ds.attrib["Calibration_stage"] == "Prior"
            @test ds.attrib["Cal_Val"] == "Calibration"
            @test ds.attrib["Conventions"] == "COARDS"
        end
        @testset "scalar + coordinate lat/lon" begin
            @test haskey(ds, "lat") && haskey(ds, "lon")
            @test haskey(ds, "latitude") && haskey(ds, "longitude")
            @test ds["latitude"].attrib["units"] == "degrees_north"
            @test ds["longitude"].attrib["units"] == "degrees_east"
        end
    end
end
