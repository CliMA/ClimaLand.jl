import ClimaTimeSteppers as CTS
using Test

problem(test_case, tab::CTS.IMEXARKAlgorithmName) = test_case.split_prob
problem(test_case, tab) = test_case.prob

algorithm(tab::CTS.IMEXARKAlgorithmName) = CTS.IMEXAlgorithm(tab, NewtonsMethod(; max_iters = 2))
algorithm(tab::CTS.RosenbrockAlgorithmName) = CTS.RosenbrockAlgorithm(ClimaTimeSteppers.tableau(tab))
algorithm(tab) = tab
