#= probably temporary script  
include("integrated/ozark/ozark.jl")
=#

using DataFrames
variables_name = ("DateTime", "RECO", "TA", "VPD", "PA", "P", "WS", "LW_IN", "SW_IN", "CO2_F", "CO2", "SWC_F",
                  "SWC", "TS_F", "TS", "GPP", "LE", "H_CORR", "H", "H_F", "G", "G_F", "LW_OUT", "SW_OUT")

variables = [vec(mat) for mat in (LOCAL_DATETIME, RECO, TA, VPD, PA, P, WS, LW_IN, SW_IN, CO2_F, CO2,
             SWC_F, SWC, TS_F, TS, GPP, LE, H_CORR, H, H_F, G, G_F, LW_OUT, SW_OUT)]

inputs = DataFrame([variables_name[i] => variables[i] for i in 1:length(variables_name)])


