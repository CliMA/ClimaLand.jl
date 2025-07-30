echo "-------------------------------- P-model at US-Var --------------------------------"
julia --project=.buildkite experiments/integrated/fluxnet/pmodel_integration_full.jl US-Var pmodel testsms1
echo "-------------------------------- Farquhar at US-Var --------------------------------"
julia --project=.buildkite experiments/integrated/fluxnet/pmodel_integration_full.jl US-Var farquhar testsms1

echo "-------------------------------- P-model at US-MOz --------------------------------"
julia --project=.buildkite experiments/integrated/fluxnet/pmodel_integration_full.jl US-MOz pmodel testsms1
echo "-------------------------------- Farquhar at US-MOz --------------------------------"
julia --project=.buildkite experiments/integrated/fluxnet/pmodel_integration_full.jl US-MOz farquhar testsms1

echo "-------------------------------- P-model at US-Ha1 --------------------------------"
julia --project=.buildkite experiments/integrated/fluxnet/pmodel_integration_full.jl US-Ha1 pmodel testsms1
echo "-------------------------------- Farquhar at US-Ha1 --------------------------------"
julia --project=.buildkite experiments/integrated/fluxnet/pmodel_integration_full.jl US-Ha1 farquhar testsms1

echo "-------------------------------- P-model at US-NR1 --------------------------------"
julia --project=.buildkite experiments/integrated/fluxnet/pmodel_integration_full.jl US-NR1 pmodel testsms1
echo "-------------------------------- Farquhar at US-NR1 --------------------------------"
julia --project=.buildkite experiments/integrated/fluxnet/pmodel_integration_full.jl US-NR1 farquhar testsms1
