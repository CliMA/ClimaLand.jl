"""
    FluxnetSimulations.make_set_fluxnet_initial_conditions(
        site_ID,
        start_date,
        hour_offset_from_UTC,
        model,
    )

Creates and returns a function `set_ic!(Y,p,t,model)` which
updates `Y` in place with an estimated set of initial conditions
based on the fluxnet observations at `site_ID` at the `start_date` in UTC,
and the type of the `model`.
In order to convert between local time and UTC, the hour offset from
UTC is required.
"""
function FluxnetSimulations.make_set_fluxnet_initial_conditions(
    site_ID,
    start_date,
    hour_offset_from_UTC,
    model,
)
    set_ic!(Y, p, t, model) =
        set_fluxnet_ic!(Y, site_ID, start_date, hour_offset_from_UTC, model)
    return set_ic!
end


"""
     set_fluxnet_ic!(
        Y,
        site_ID,
        start_date,
        hour_offset_from_UTC,
        model::ClimaLand.AbstractLandModel,
    )

Sets the initial conditions of `Y` using observations from the site `site_ID`, if available,
using the observations closest to the start_date (in UTC). Since the data from Fluxnet sites
is provided in local time, we require the offset from UTC in hours `hour_offset_from_UTC`.
The `model` indicates which how to update it `Y` from these observations,
via different methods of `set_fluxnet_ic!`.
"""
function set_fluxnet_ic!(
    Y,
    site_ID,
    start_date,
    hour_offset_from_UTC,
    model::ClimaLand.AbstractLandModel,
)
    fluxnet_csv_path = ClimaLand.Artifacts.experiment_fluxnet_data_path(site_ID)

    # Read the data and get the column name map
    (data, columns) = readdlm(fluxnet_csv_path, ','; header = true)

    # Get the UTC datetimes of the data
    varnames = ("TIMESTAMP_START",)
    column_name_map =
        get_column_name_map(varnames, columns; error_on_missing = true)
    UTC_datetimes = get_UTC_datetimes(
        hour_offset_from_UTC,
        data,
        column_name_map;
        timestamp_name = "TIMESTAMP_START",
    )
    Δ_date = UTC_datetimes .- start_date
    for component in ClimaLand.land_components(model)
        set_fluxnet_ic!(Y, data, columns, Δ_date, getproperty(model, component))
    end
end

"""
     set_fluxnet_ic!(Y, data, columns, Δ_date, model::ClimaLand.Soil.EnergyHydrology)

Sets the values of Y.soil in place with:
- \\vartheta_l: a profile `θ_deep + (θ_sfc - θ_deep) exp(z / z_θ)` (`z ≤ 0`, `z_θ = 0.5 m`)
  between the observed shallow water content `θ_sfc` at the observation date closest to the
  start date and a deep value `θ_deep`, the long-term mean of the deepest available soil
  moisture record (`SWC_F_MDS_2`, else `SWC_F_MDS_1`). The seasonal moisture signal is
  confined to roughly the top meter of soil, below which the water content is close to its
  long-term mean; starting the whole column at the shallow value of one date puts, at a
  semi-arid site, meters of water into a column that only loses water by evapotranspiration.
  Where the site has no soil moisture record at all, `θ_deep` is estimated from the
  climate ([`climatological_soil_moisture`](@ref)) and `θ_sfc = θ_deep`. All values are
  bounded between the permanent wilting point (ψ = -150 m) and 95% of the effective
  saturation range above the residual water content.
- θ_i: no ice (θ_i = 0)
- \\rho e_int: an internal energy computed using the above θ_l, θ_i, and a temperature
  profile `T_deep + (T_sfc - T_deep) exp(z / z_T)` (`z_T = 2 m`, the annual damping depth)
  between the shallow soil temperature at the observation date closest to the start date
  (the air temperature if unavailable) and the record mean of the same column.

Here, `Y` is the prognostic field vector, `data` is the raw data for the site read from
a CSV file, `columns` is the list of column names,
`Δ_date` is the vector of date differences between the observations (in UTC) and the
start date (in UTC), and `model` indicates which part of `Y` we are updating, and how to update it,
via different methods of `set_fluxnet_ic!`.
"""
function set_fluxnet_ic!(
    Y,
    data,
    columns,
    Δ_date,
    model::ClimaLand.Soil.EnergyHydrology;
    val = -9999,
)
    FT = eltype(Y.soil.ρe_int)
    (; θ_r, ν, hydrology_cm) = model.parameters
    column(name) = fluxnet_column(data, columns, name; val)
    record_mean(v) = sum(x for x in v if x != val) / count(!=(val), v)

    # Soil moisture: shallow sensor at the start date, deepest sensor mean at depth
    swc_1 = column("SWC_F_MDS_1")
    swc_2 = column("SWC_F_MDS_2")
    ts_1 = column("TS_F_MDS_1")
    if !isnothing(swc_1) && !isnothing(ts_1)
        # Frozen records are masked with the missing-value marker
        unfrozen_swc = ifelse.(ts_1 .> 0, swc_1, oftype(first(swc_1), val))
        if any(x -> !var_missing(x; val), unfrozen_swc)
            swc_1 = unfrozen_swc
        end
    end
    swc_sfc = isnothing(swc_1) ? swc_2 : swc_1
    swc_deep = isnothing(swc_2) ? swc_1 : swc_2
    # Bounds of the hydraulics: permanent wilting point and 95% of saturation
    θ_wilt = @. θ_r +
       (ν - θ_r) *
       ClimaLand.Soil.inverse_matric_potential(hydrology_cm, FT(-150))
    θ_max = @. θ_r + (ν - θ_r) * FT(0.95)
    if isnothing(swc_sfc)
        θ_fc = @. θ_r +
           (ν - θ_r) *
           ClimaLand.Soil.inverse_matric_potential(hydrology_cm, FT(-3.3))
        θ_deep = climatological_soil_moisture(data, columns, θ_wilt, θ_fc; val)
        θ_sfc = θ_deep
    else
        θ_sfc = FT(
            get_data_at_start_date(
                swc_sfc,
                Δ_date;
                preprocess_func = x -> x / 100,
                val,
                varname = "SWC_F_MDS",
            ),
        )
        θ_deep = FT(record_mean(swc_deep) / 100)
    end
    z = model.domain.fields.z
    z_θ = FT(0.5)
    # The retention curve is defined for ϑ_l > θ_r only. Where the observed
    # water content is below the residual water content of the soil
    # parameters, the soil is initialized at the water content of the
    # permanent wilting point (ψ = -150 m), the driest state the hydraulics
    # represent: at ϑ_l ≤ θ_r the pressure head is unbounded while its
    # derivative vanishes, and the first wetting of the surface layer then
    # drives an unbounded flux.
    @. Y.soil.ϑ_l =
        clamp(θ_deep + (θ_sfc - θ_deep) * exp(z / z_θ), θ_wilt, θ_max)
    Y.soil.θ_i .= 0

    # Soil temperature: shallow sensor at the start date, record mean at depth
    T_col = isnothing(ts_1) ? column("TA_F") : ts_1
    T_sfc = FT(
        get_data_at_start_date(
            T_col,
            Δ_date;
            preprocess_func = x -> x + 273.15,
            val,
            varname = isnothing(ts_1) ? "TA_F" : "TS_F_MDS_1",
        ),
    )
    T_deep = FT(record_mean(T_col) + 273.15)
    z_T = FT(2)
    T_soil_0 = @. T_deep + (T_sfc - T_deep) * exp(z / z_T)

    ρc_s = ClimaLand.Soil.volumetric_heat_capacity.(
        Y.soil.ϑ_l,
        Y.soil.θ_i,
        model.parameters.ρc_ds,
        model.parameters.earth_param_set,
    )
    Y.soil.ρe_int = ClimaLand.Soil.volumetric_internal_energy.(
        Y.soil.θ_i,
        ρc_s,
        T_soil_0,
        model.parameters.earth_param_set,
    )
end

"""
    fluxnet_column(data, columns, name; val = -9999)

Return the column `name` of the FLUXNET `data` matrix, or `nothing` if the column is
absent or was never observed at the site (all entries equal to `val`).
"""
function fluxnet_column(data, columns, name; val = -9999)
    idx = findfirst(columns[:] .== name)
    (isnothing(idx) || all_missing(data[:, idx]; val)) && return nothing
    return data[:, idx]
end

"""
    climatological_soil_moisture(data, columns, θ_wilt, θ_fc; val = -9999)

Estimate the long-term soil water content at a site without a soil moisture record from
its climate, as `θ_wilt + (θ_fc - θ_wilt) min(1, P/PET)`: a wet climate (`P ≥ PET`) holds
the soil near field capacity `θ_fc`, a dry one near the wilting point `θ_wilt`. `P` is the
record-mean precipitation rate (`P_F`, mm per half hour) and `PET` the Priestley-Taylor
potential evaporation `1.26 Δ/(Δ+γ) R_n/λ` with `Δ/(Δ+γ) = 0.65` (about 15 °C) and the
record-mean net radiation (`NETRAD`, or `0.55 SW_IN_F` where net radiation was not
measured). `θ_wilt` and `θ_fc` and the return value are fields on the soil domain.
"""
function climatological_soil_moisture(data, columns, θ_wilt, θ_fc; val = -9999)
    FT = eltype(θ_fc)
    valid_mean(v) = sum(x for x in v if x != val) / count(!=(val), v)
    P = fluxnet_column(data, columns, "P_F"; val)
    isnothing(P) && error(
        "FLUXNET site has neither soil moisture nor precipitation records, so the \
         soil initial condition cannot be estimated.",
    )
    precip = valid_mean(P) / 1800 # mm per half hour -> mm/s
    Rn = fluxnet_column(data, columns, "NETRAD"; val)
    SW = fluxnet_column(data, columns, "SW_IN_F"; val)
    net_rad = isnothing(Rn) ? 0.55 * valid_mean(SW) : valid_mean(Rn) # W/m²
    λ = 2.5e6 # J/kg; 1 mm of water is 1 kg/m²
    pet = max(1.26 * 0.65 * net_rad / λ, eps(Float64)) # mm/s
    aridity = FT(min(1, precip / pet))
    @info "Soil moisture initial condition estimated from climate" precip *
                                                                   86400 pet *
                                                                         86400 aridity
    return @. θ_wilt + (θ_fc - θ_wilt) * aridity
end

"""
    set_fluxnet_ic!(Y, data, columns, Δ_date, model::ClimaLand.Canopy.CanopyModel)

Sets Y.canopy.energy.T to the air temperature at the observation date closest to the start
date of the model; sets the potential in the stem and leaf to -0.1 and -0.2 MPa, respectively,
and the computes the resulting water content Y.canopy.hydraulics.ϑ_l using the retention curve
of the plant. If the biomass model carries prognostic state (the optimal-LAI model), that is
set from its climatology.
"""
function set_fluxnet_ic!(
    Y,
    data,
    columns,
    Δ_date,
    model::ClimaLand.Canopy.CanopyModel{FT};
    val = -9999,
) where {FT}
    if model.energy isa Canopy.BigLeafEnergyModel
        idx = findfirst(columns[:] .== "TA_F")
        T_air_0 = get_data_at_start_date(
            data[:, idx],
            Δ_date;
            preprocess_func = x -> x + 273.15,
            val,
            varname = "TA_F",
        )
        Y.canopy.energy.T .= T_air_0
    end

    ψ_leaf_0 = FT(-2e5 / 9800)
    hydraulics = model.hydraulics
    S_l_ini = ClimaLand.Canopy.inverse_water_retention_curve.(
        hydraulics.parameters.retention_model,
        ψ_leaf_0,
        hydraulics.parameters.ν,
        hydraulics.parameters.S_s,
    )
    Y.canopy.hydraulics.ϑ_l .= ClimaLand.Canopy.augmented_liquid_fraction.(
        hydraulics.parameters.ν,
        S_l_ini,
    )
    ClimaLand.Simulations.set_canopy_component_initial_conditions!(
        Y,
        nothing,
        model.biomass,
        model,
    )
    if model.photosynthesis isa Canopy.PModel
        idx = findfirst(columns[:] .== "VPD_F")
        VPD_0 = get_data_at_start_date(
            data[:, idx],
            Δ_date;
            preprocess_func = x -> x * 100, # hPa to Pa
            val,
        )
        idx = findfirst(columns[:] .== "PA_F")
        P_air_0 = get_data_at_start_date(
            data[:, idx],
            Δ_date;
            preprocess_func = x -> x * 1000,# kPa to Pa
            val,
        )
        idx = findfirst(columns[:] .== "CO2_F_MDS")
        c_co2_0 = get_data_at_start_date(
            data[:, idx],
            Δ_date;
            preprocess_func = x -> x * 1e-6, # convert from μmol/mol to mol/mol
            val,
        )
        idx = findfirst(columns[:] .== "TA_F")
        T_air_0 = get_data_at_start_date(
            data[:, idx],
            Δ_date;
            preprocess_func = x -> x + 273.15,# C to K
            val,
        )
        βm = FT(1)
        APAR_canopy_moles = FT(1e-3) # mol/m^2/s, from 1000 μmol/m^2/s
        @. Y.canopy.photosynthesis.acclimated = compute_optimal_capacities(
            model.photosynthesis.parameters,
            model.photosynthesis.constants,
            FT(T_air_0),
            FT(P_air_0),
            FT(VPD_0),
            FT(c_co2_0),
            βm,
            APAR_canopy_moles,
        )
    end

end

"""
    set_fluxnet_ic!(Y, data, columns, Δ_date, model::ClimaLand.Snow.SnowModel)

Sets Y.snow.S, Y.snow.S_l, and Y.snow.U in place to be zero at the start of the simulation
(no snow).

Note that the Snow NeuralDensity model has additional prognostic variables which also must be set
to zero; another method may work well for that case.
"""
function set_fluxnet_ic!(
    Y,
    data,
    columns,
    Δ_date,
    model::ClimaLand.Snow.SnowModel,
)
    Y.snow.S .= 0.0
    Y.snow.S_l .= 0.0
    Y.snow.U .= 0.0
end

"""
     set_fluxnet_ic!(Y, data, columns, Δ_date, model::ClimaLand.Soil.Biogeochemistry.SoilCO2Model)

Sets Y.soilco2.CO2, Y.soilco2.O2, and Y.soilco2.SOC in place with initial values.

The CO2 and O2 prognostic variables are stored as total mass per unit *bulk soil
volume* (storage = θ_eff · gas-phase concentration, where θ_eff = θ_a + β·θ_l
accounts for both the air-filled pore space and the gas dissolved in pore water):
- CO2: total CO₂ as carbon mass per bulk soil volume (kg C m⁻³ soil)
- O2: total O₂ (gas-phase in air-filled pores + dissolved in pore water) per bulk
  soil volume (kg O₂ m⁻³ soil)
- SOC: soil organic carbon concentration (kg C/m³)

The CO2/O2 values are representative atmospheric-equilibrium initial guesses
(assuming θ_eff ≈ 0.3 at standard T, P); the near-surface values are rapidly
replaced by the diffusive boundary condition.
"""
function set_fluxnet_ic!(
    Y,
    data,
    columns,
    Δ_date,
    model::ClimaLand.Soil.Biogeochemistry.SoilCO2Model,
)
    # Representative bulk densities (kg m⁻³ soil) ≈ θ_eff · gas-phase density.
    Y.soilco2.CO2 .= 6e-5   # ≈ θ_eff · c_atm · P·M_C/(R·T), kg C m⁻³ soil
    Y.soilco2.O2 .= 0.08    # ≈ θ_eff · 0.21 · P·M_O2/(R·T), kg O₂ m⁻³ soil
    Y.soilco2.SOC .= 5.0    # Default SOC concentration (kg C/m³)
end
