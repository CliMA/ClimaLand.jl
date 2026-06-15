#!/usr/bin/env bash
set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
STATE_DIR="$REPO_DIR/.automation/callmip_watchdog"
LOG_FILE="$STATE_DIR/watchdog.log"
HANDLED_FILE="$STATE_DIR/handled_failed_jobs.txt"
ALERT_EMAIL="${WATCHDOG_ALERT_EMAIL:-renato.braghiere@gmail.com}"
WEBHOOK_URL="${WATCHDOG_WEBHOOK_URL:-${SLACK_WEBHOOK_URL:-}}"
HOST_ID="$(hostname -f 2>/dev/null || hostname)"
mkdir -p "$STATE_DIR"
touch "$LOG_FILE" "$HANDLED_FILE"

log() {
  local msg="$1"
  printf '%s %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$msg" >>"$LOG_FILE"
}

json_escape() {
  local s="$1"
  s=${s//\\/\\\\}
  s=${s//\"/\\\"}
  s=${s//$'\n'/\\n}
  printf '%s' "$s"
}

send_email_alert() {
  local subject="$1"
  local body="$2"

  if [[ -z "$ALERT_EMAIL" ]]; then
    return 0
  fi

  if command -v mail >/dev/null 2>&1; then
    printf '%s\n' "$body" | mail -s "$subject" "$ALERT_EMAIL" || log "alert email failed via mail"
  elif command -v mailx >/dev/null 2>&1; then
    printf '%s\n' "$body" | mailx -s "$subject" "$ALERT_EMAIL" || log "alert email failed via mailx"
  else
    log "alert email skipped: mail/mailx not available"
  fi
}

send_webhook_alert() {
  local text="$1"

  if [[ -z "$WEBHOOK_URL" ]]; then
    return 0
  fi

  local payload
  payload=$(printf '{"text":"%s"}' "$(json_escape "$text")")
  curl -fsS -X POST -H 'Content-Type: application/json' --data "$payload" "$WEBHOOK_URL" >/dev/null || \
    log "alert webhook failed"
}

send_alert() {
  local level="$1"
  local message="$2"
  local details="${3:-}"
  local subject="[callmip-watchdog][$level] $HOST_ID"
  local body

  body="$(cat <<EOF
$subject
time: $(date '+%Y-%m-%d %H:%M:%S %Z')
repo: $REPO_DIR
message: $message

$details
EOF
)"

  send_email_alert "$subject" "$body"
  send_webhook_alert "$body"
}

cd "$REPO_DIR"

# Prevent overlapping runs.
exec 9>"$STATE_DIR/lock"
if ! flock -n 9; then
  log "another watchdog run is active; skipping"
  exit 0
fi

apply_known_fixes() {
  local changed=0

  # Fix known stale canopy temperature field in autotrophic respiration.
  if grep -q "p.canopy.energy.T" src/standalone/Vegetation/autotrophic_respiration.jl; then
    sed -i 's/p\.canopy\.energy\.T/p.drivers.T/g' src/standalone/Vegetation/autotrophic_respiration.jl
    changed=1
    log "applied fix: autotrophic respiration now uses p.drivers.T"
  fi

  # Fix known ParamDict haskey pattern in autotrophic respiration.
  if grep -q 'haskey(toml_dict, "autotrophic_respiration_Q10")' src/standalone/Vegetation/autotrophic_respiration.jl; then
    perl -0777 -i -pe 's/Q10 = haskey\(toml_dict, "autotrophic_respiration_Q10"\) \?\n\s*toml_dict\["autotrophic_respiration_Q10"\] :\n\s*CP\.float_type\(toml_dict\)\(2\.0\),/Q10 = try\n        toml_dict["autotrophic_respiration_Q10"]\n    catch err\n        err isa KeyError || rethrow(err)\n        CP.float_type(toml_dict)(2.0)\n    end,/s' src/standalone/Vegetation/autotrophic_respiration.jl
    changed=1
    log "applied fix: autotrophic respiration ParamDict optional key handling"
  fi

  # Fix known ParamDict haskey pattern in soil biogeochemistry.
  if grep -q 'haskey(toml_dict, "soilCO2_reference_rate")' src/standalone/Soil/Biogeochemistry/Biogeochemistry.jl; then
    perl -0777 -i -pe 's/\# Accept either legacy α_sx or centered Arrhenius soilCO2_reference_rate\.\n\s*α_sx = if haskey\(toml_dict, "soilCO2_reference_rate"\)\n\s*R = FT\(LP\.gas_constant\(earth_param_set\)\)\n\s*T_ref_sx = FT\(288\.15\)\n\s*V_ref_sx = FT\(toml_dict\["soilCO2_reference_rate"\]\)\n\s*V_ref_sx \/ exp\(-parameters\.Ea_sx \/ \(R \* T_ref_sx\)\)\n\s*elseif haskey\(toml_dict, "soilCO2_pre_exponential_factor"\)\n\s*FT\(toml_dict\["soilCO2_pre_exponential_factor"\]\)\n\s*else\n\s*error\("Missing soil CO2 rate parameter: expected soilCO2_reference_rate or soilCO2_pre_exponential_factor"\)\n\s*end/\# Accept either legacy α_sx or centered Arrhenius soilCO2_reference_rate.\n    α_sx = begin\n        reference_rate = try\n            toml_dict["soilCO2_reference_rate"]\n        catch err\n            err isa KeyError || rethrow(err)\n            nothing\n        end\n\n        if reference_rate !== nothing\n            R = FT(LP.gas_constant(earth_param_set))\n            T_ref_sx = FT(288.15)\n            V_ref_sx = FT(reference_rate)\n            V_ref_sx \/ exp(-parameters.Ea_sx \/ (R * T_ref_sx))\n        else\n            pre_exponential = try\n                toml_dict["soilCO2_pre_exponential_factor"]\n            catch err\n                err isa KeyError || rethrow(err)\n                nothing\n            end\n\n            if pre_exponential !== nothing\n                FT(pre_exponential)\n            else\n                error("Missing soil CO2 rate parameter: expected soilCO2_reference_rate or soilCO2_pre_exponential_factor")\n            end\n        end\n    end/s' src/standalone/Soil/Biogeochemistry/Biogeochemistry.jl
    changed=1
    log "applied fix: soil biogeochemistry ParamDict optional key handling"
  fi

  echo "$changed"
}

# If calibration is still active, do nothing.
if squeue -u "$USER" -h -o "%j|%T" | awk -F'|' '$1 ~ /dk_sor_calibration/ && ($2 == "RUNNING" || $2 == "PENDING") { found=1 } END { exit(found ? 0 : 1) }'; then
  log "calibration is active; no action"
  exit 0
fi

lookup_start="$(date -d '2 days ago' +%F)"
latest_failed_job="$(sacct -u "$USER" --starttime "$lookup_start" --format=JobIDRaw,JobName,State -n -X | awk '$2 ~ /dk_sor_ca/ && $3 ~ /FAILED/ { print $1 }' | tail -n 1)"
if [[ -z "$latest_failed_job" ]]; then
  log "no failed calibration job found since $lookup_start"
  exit 0
fi

if grep -qx "$latest_failed_job" "$HANDLED_FILE"; then
  log "failed job $latest_failed_job already handled"
  exit 0
fi

err_file="$REPO_DIR/dk_sor_calibration_${latest_failed_job}.err"
out_file="$REPO_DIR/dk_sor_calibration_${latest_failed_job}.out"
if [[ -f "$err_file" ]]; then
  log "handling failed calibration job $latest_failed_job"
  tail -n 120 "$err_file" >"$STATE_DIR/last_error_${latest_failed_job}.log" || true
  err_excerpt="$(tail -n 40 "$STATE_DIR/last_error_${latest_failed_job}.log" 2>/dev/null || true)"
  send_alert "WARN" "Detected failed calibration job $latest_failed_job" "Error tail file: $STATE_DIR/last_error_${latest_failed_job}.log

Recent error excerpt:
$err_excerpt"
else
  log "failed job $latest_failed_job has no local err file"
  send_alert "WARN" "Failed calibration job $latest_failed_job has no local err file"
fi

changed="$(apply_known_fixes)"

if [[ "$changed" -eq 1 ]]; then
  log "known fixes applied; running readiness gate"
else
  log "no known fixes needed; running readiness gate"
fi

if ./test_callmip_pipeline.sh >"$STATE_DIR/readiness_${latest_failed_job}.log" 2>&1; then
  log "readiness gate passed; cancelling dead dependency tail"
  dep_ids="$(squeue -u "$USER" -h -o "%i|%j|%T|%R" | awk -F'|' '$2 ~ /(dk_sor_emulate_sample|callmip_sims|callmip_postprocess)/ && $3 == "PENDING" && $4 ~ /Dependency/ { print $1 }')"
  if [[ -n "$dep_ids" ]]; then
    # shellcheck disable=SC2086
    scancel $dep_ids || true
    log "cancelled dependency-pending jobs: $dep_ids"
  fi

  log "submitting fresh pipeline chain"
  bash submit_callmip_pipeline.sh >"$STATE_DIR/resubmit_${latest_failed_job}.log" 2>&1
  log "resubmission complete"
  send_alert "INFO" "Resubmitted pipeline after failed job $latest_failed_job" "Submission log: $STATE_DIR/resubmit_${latest_failed_job}.log"
else
  log "readiness gate failed; manual intervention required. See $STATE_DIR/readiness_${latest_failed_job}.log"
  readiness_excerpt="$(tail -n 60 "$STATE_DIR/readiness_${latest_failed_job}.log" 2>/dev/null || true)"
  send_alert "ERROR" "Readiness gate failed after failed calibration job $latest_failed_job" "Readiness log: $STATE_DIR/readiness_${latest_failed_job}.log

Recent readiness excerpt:
$readiness_excerpt"
fi

echo "$latest_failed_job" >>"$HANDLED_FILE"
log "marked failed job $latest_failed_job as handled"
