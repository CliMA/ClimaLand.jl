#= run those lines for testing / development

ARGS = ["US-MOz"]
include("integrated/fluxnet/setup.jl")

=#

using DataFrames
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
        LOCAL_DATETIME,
        drivers.RECO.values,
        drivers.TA.values,
        drivers.VPD.values,
        drivers.PA.values,
        drivers.P.values,
        drivers.WS.values,
        drivers.LW_IN.values,
        drivers.SW_IN.values,
        drivers.CO2.values,
        drivers.SWC.values,
        drivers.TS.values,
        drivers.GPP.values,
        drivers.LE.values,
        drivers.H.values,
        drivers.G.values,
        drivers.SW_OUT.values,
    )
]

inputs = DataFrame([
    variables_name[i] => variables[i] for i in 1:length(variables_name)
])

# make a Unitful dataframe with SI units and one with commonly used units

using Unitful: R, L, mol, K, kJ, °C, m, g, cm, hr, mg, s, μmol, Pa, W, mm
using UnitfulMoles: molC
using Unitful, UnitfulMoles
@compound CO₂

columns = variables_name[2:end]
units = [molCO₂*m^-2*s^-1, K, Pa, Pa, m*s^-1, m*s^-1, W*m^-2, W*m^-2, molCO₂, m^3*m^-3, K, molCO₂*m^-2*s^-1, W*m^-2, W*m^-2, W*m^-2, W*m^-2]
inputs_SI = copy(inputs)
foreach((col, unit) -> inputs_SI[!, col] .= first([inputs_SI[!, col]]unit), columns, units)

units_to = [μmolCO₂*m^-2*s^-1, °C, kPa, Pa, mm*s^-1, m*s^-1, W*m^-2, W*m^-2, μmolCO₂, m^3*m^-3, °C, μmolCO₂*m^-2*s^-1, W*m^-2, W*m^-2, W*m^-2, W*m^-2]
inputs_commonly_used = copy(inputs_SI)
foreach((col, unit_to) -> inputs_commonly_used[!, col] = uconvert.(unit_to, inputs_SI[!, col]), columns, units_to)


