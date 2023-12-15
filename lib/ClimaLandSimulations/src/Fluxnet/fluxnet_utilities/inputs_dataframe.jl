export make_inputs_df

# UnitfulMoles compound
@compound CO2

"""
    make_inputs_df(site_ID)

Returns DataFrames of inputs, one unitless, one in SI units, one in commonly used units. 
"""
function make_inputs_df(
    site_ID;
    FT = Float64,
    context = ClimaComms.SingletonCommsContext(),
    setup = make_setup(site_ID),
    domain = make_domain(setup, FT),
    config = make_config(site_ID),
    params = make_parameters(site_ID),
    drivers = make_drivers(site_ID, setup, config, params, context),
    timestepper = make_timestepper(setup),
)

    variables_name = (
        "DateTime",
        "RECO",
        "TA",
        "VPD",
        "PA",
        "P",
        "WS",
        "LW_IN",
        "SW_IN",
        "CO2",
        "SWC",
        "TS",
        "GPP",
        "LE",
        "H",
        "G",
        "SW_OUT",
    ) # I imagine this may not work for all fluxnet sites...

    variables = [
        vec(mat) for mat in (
            drivers.LOCAL_DATETIME,
            drivers.drivers.RECO.values,
            drivers.drivers.TA.values,
            drivers.drivers.VPD.values,
            drivers.drivers.PA.values,
            drivers.drivers.P.values,
            drivers.drivers.WS.values,
            drivers.drivers.LW_IN.values,
            drivers.drivers.SW_IN.values,
            drivers.drivers.CO2.values,
            drivers.drivers.SWC.values,
            drivers.drivers.TS.values,
            drivers.drivers.GPP.values,
            drivers.drivers.LE.values,
            drivers.drivers.H.values,
            drivers.drivers.G.values,
            drivers.drivers.SW_OUT.values,
        )
    ]

    inputs = DataFrame([
        variables_name[i] => variables[i] for i in 1:length(variables_name)
    ])

    # make a Unitful dataframe with SI units and one with commonly used units
    columns = variables_name[2:end]
    units = [
        molCO₂ * m^-2 * s^-1,
        K,
        Pa,
        Pa,
        m * s^-1,
        m * s^-1,
        W * m^-2,
        W * m^-2,
        molCO₂,
        m^3 * m^-3,
        K,
        molCO₂ * m^-2 * s^-1,
        W * m^-2,
        W * m^-2,
        W * m^-2,
        W * m^-2,
    ]
    inputs_SI = copy(inputs)
    foreach(
        (col, unit) -> inputs_SI[!, col] .= first([inputs_SI[!, col]]unit),
        columns,
        units,
    )

    units_to = [
        μmolCO₂ * m^-2 * s^-1,
        °C,
        kPa,
        Pa,
        mm * s^-1,
        m * s^-1,
        W * m^-2,
        W * m^-2,
        μmolCO₂,
        m^3 * m^-3,
        °C,
        μmolCO₂ * m^-2 * s^-1,
        W * m^-2,
        W * m^-2,
        W * m^-2,
        W * m^-2,
    ]
    inputs_commonly_used = copy(inputs_SI)
    foreach(
        (col, unit_to) ->
            inputs_commonly_used[!, col] =
                uconvert.(unit_to, inputs_SI[!, col]),
        columns,
        units_to,
    )

    return (inputs, inputs_SI, inputs_commonly_used)

end
