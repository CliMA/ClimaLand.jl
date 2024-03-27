Y, p, cds = initialize(land)
#Initial conditions
Y.soil.ϑ_l = drivers.SWC.values[1 + Int(round(t0 / 1800))] # Get soil water content at t0
# recalling that the data is in intervals of 1800 seconds. Both the data
# and simulation are reference to 2005-01-01-00 (LOCAL)
# or 2005-01-01-06 (UTC)
Y.soil.θ_i = FT(0.0)
T_0 = FT(drivers.TS.values[1 + Int(round(t0 / 1800))]) # Get soil temperature at t0
ρc_s =
    volumetric_heat_capacity.(Y.soil.ϑ_l, Y.soil.θ_i, Ref(land.soil.parameters))
Y.soil.ρe_int =
    volumetric_internal_energy.(
        Y.soil.θ_i,
        ρc_s,
        T_0,
        Ref(land.soil.parameters),
    )

Y.soilco2.C .= FT(0.000412) # set to atmospheric co2, mol co2 per mol air

ψ_stem_0 = FT(-1e5 / 9800)
ψ_leaf_0 = FT(-2e5 / 9800)

S_l_ini =
    inverse_water_retention_curve.(
        retention_model,
        [ψ_stem_0, ψ_leaf_0],
        plant_ν,
        plant_S_s,
    )

for i in 1:2
    Y.canopy.hydraulics.ϑ_l.:($i) .=
        augmented_liquid_fraction.(plant_ν, S_l_ini[i])
end

Y.canopy.energy.T = drivers.TA.values[1 + Int(round(t0 / 1800))] # Get atmos temperature at t0
set_initial_cache!(p, Y, t0)
