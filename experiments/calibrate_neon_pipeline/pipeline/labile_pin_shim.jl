"""
labile_pin_shim.jl — defines the METHODS that make a SINGLE model interface serve
both "labile calibrated" and "labile off" modes.

The `NeonPipelineInterface` STRUCT itself is defined at top level in
Calibration.jl (so its constructor isn't world-age-too-new when the main process
builds it). This file only adds the forward_model / observation_map methods, and
is loaded — together with the original `model_interface_wLabile.jl` — on the
workers (and the main process) via the broadcast in run_calibration. By then
`NeonLabileModelInterface`, the original `forward_model`/`observation_map`, and
the worker globals (SITE_ID, OUTPUT_DIR, OBS_FILEPATH, DT, Caldepthnum,
LABILE_ON) are all defined.

`model_interface_wLabile.jl` reads `labile_depth_scale` directly from each
per-member parameter TOML (it is not a ClimaParams key). When the parameter is
NOT calibrated, ClimaCalibrate never writes that key, so the read would throw
`KeyError`. `forward_model` below first ensures the per-member TOML contains
`labile_depth_scale` (injecting 0.0 when absent — `exp(0·z)=1`, the plain
model), then delegates to `forward_model(::NeonLabileModelInterface, …)`.
`observation_map` simply forwards to the original.
"""

import ClimaCalibrate
import TOML as _PIPELINE_TOML

# Dispatch tag. Defined here (loaded at top level by Calibration.jl, and broadcast
# to workers) so its constructor + methods live in an old-enough world age.
struct NeonPipelineInterface <: ClimaCalibrate.AbstractModelInterface end

# Ensure the per-member parameter TOML has a labile_depth_scale entry. When the
# parameter is calibrated, EKP already wrote it (no-op); otherwise inject k=0.
function _ensure_labile_entry!(iteration, member)
    path = ClimaCalibrate.parameter_path(OUTPUT_DIR, iteration, member)
    data = _PIPELINE_TOML.parsefile(path)
    if !haskey(data, "labile_depth_scale")
        data["labile_depth_scale"] =
            Dict("value" => 0.0, "type" => "float", "used_in" => ["Land"])
        open(path, "w") do io
            _PIPELINE_TOML.print(io, data)
        end
    end
    return nothing
end

function ClimaCalibrate.forward_model(
    ::NeonPipelineInterface, iteration, member,
)
    _ensure_labile_entry!(iteration, member)
    # Delegate to the original wLabile implementation.
    return ClimaCalibrate.forward_model(
        NeonLabileModelInterface(), iteration, member,
    )
end

function ClimaCalibrate.observation_map(::NeonPipelineInterface, iteration)
    return ClimaCalibrate.observation_map(NeonLabileModelInterface(), iteration)
end
