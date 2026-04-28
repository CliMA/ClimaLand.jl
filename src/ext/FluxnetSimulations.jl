module FluxnetSimulations

# FLUXNET columns that `prescribed_forcing_fluxnet` builds TimeVaryingInputs from.
# Exposed here (rather than only in the ext) so callers that need to keep an
# observation/simulation window in sync — e.g. site-level calibration — can pass
# this list to `get_data_dates(...; required_columns)` and get the same advanced
# start_date the forward model uses internally.
const FLUXNET_FORCING_COLUMNS = (
    "TA_F",
    "VPD_F",
    "PA_F",
    "P_F",
    "WS_F",
    "LW_IN_F",
    "SW_IN_F",
    "CO2_F_MDS",
)

function prescribed_forcing_fluxnet end

function prescribed_LAI_fluxnet end

function make_set_fluxnet_initial_conditions end

function get_data_dt end

function get_comparison_data end

function get_data_dates end

function get_maxLAI_at_site end

function get_domain_info end

function get_location end

function get_fluxtower_height end

function get_parameters end

function replace_hyphen end

# Generic site driver and FLUXNET2015 metadata helpers.
function generic_site_simulation end

function get_site_info end

function get_canopy_height end

end
