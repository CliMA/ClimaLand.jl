module FluxnetSimulations

# FLUXNET columns that `prescribed_forcing_fluxnet` builds TimeVaryingInputs from.
# Pass to `get_data_dates(...; required_columns)` so that the forcing covers the
# whole simulation.
const FLUXNET_FORCING_COLUMNS =
    ("TA_F", "VPD_F", "PA_F", "P_F", "WS_F", "LW_IN_F", "SW_IN_F", "CO2_F_MDS")

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

# FLUXNET2015 metadata helpers.
function get_site_info end

function get_canopy_height end

function get_site_igbp end

end
