SITES=(
  "US-Var"
  "US-NR1"
  "BR-Sa3"
  "AU-Emr"
  "US-MMS"
  "US-IB2"
  "AU-Rob"
  "CN-Cha"
  "US-SRM"
  "AU-Cpr"
)

for site in "${SITES[@]}"; do
    echo "-------------------------------- P-model piecewise sms $site --------------------------------"
    julia --project=.buildkite experiments/standalone/Vegetation/pmodel_integration.jl "$site" pmodel piecewise
    echo "-------------------------------- P-model no sms $site --------------------------------"
    julia --project=.buildkite experiments/standalone/Vegetation/pmodel_integration.jl "$site" pmodel no_sms
    echo "-------------------------------- Farquhar piecewise sms $site --------------------------------"
    julia --project=.buildkite experiments/standalone/Vegetation/pmodel_integration.jl "$site" farquhar piecewise
    echo "-------------------------------- Farquhar no sms $site --------------------------------"
    julia --project=.buildkite experiments/standalone/Vegetation/pmodel_integration.jl "$site" farquhar no_sms
done