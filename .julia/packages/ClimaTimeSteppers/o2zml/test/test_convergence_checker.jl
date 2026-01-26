using ClimaTimeSteppers, Test
import ClimaTimeSteppers as CTS

@testset "ConvergenceChecker" begin
    val_func(iter) = [60.0, -80.0]
    err_func1(iter) = [6, -8] .* (-1 / 10)^min(iter, 11)
    err_func2(iter) = [60, -80] .* (-1 / 10)^min(iter, 12)
    conditions = (
        MaximumError(1e-10 + eps()),
        MaximumRelativeError(1e-10 + eps()),
        MaximumErrorReduction(1e-10 + eps()),
        MinimumRateOfConvergence(0.1 + eps(), 1),
    )
    last_iters1 = (11, 9, 10, 12)
    last_iters2 = (12, 10, 10, 13)
    function test_func(checker, last_iter1, last_iter2)
        cache = CTS.allocate_cache(checker, val_func(0))
        for (err_func, last_iter) in ((err_func1, last_iter1), (err_func2, last_iter2))
            for iter in 0:(last_iter - 1)
                CTS.is_converged!(checker, cache, val_func(iter), err_func(iter), iter) && return false
            end
            CTS.is_converged!(checker, cache, val_func(last_iter), err_func(last_iter), last_iter) || return false
        end
        return true
    end
    for (condition, last_iter1, last_iter2) in zip(conditions, last_iters1, last_iters2)
        norm_checker = ConvergenceChecker(; norm_condition = condition)
        component_checker = ConvergenceChecker(; component_condition = condition)
        for checker in (norm_checker, component_checker)
            @test test_func(checker, last_iter1, last_iter2)
        end
    end
    for (condition_combiner, last_iter_combiner) in ((all, maximum), (any, minimum))
        condition = MultipleConditions(condition_combiner, conditions...)
        last_iter1 = last_iter_combiner(last_iters1)
        last_iter2 = last_iter_combiner(last_iters2)
        norm_checker = ConvergenceChecker(; norm_condition = condition)
        component_checker = ConvergenceChecker(; component_condition = condition)
        for checker in (norm_checker, component_checker)
            @test test_func(checker, last_iter1, last_iter2)
        end
    end
    for (condition_combiner, last_iter_combiner) in ((&, maximum), (|, minimum))
        all_condition = MultipleConditions(all, conditions...)
        any_condition = MultipleConditions(any, conditions...)
        last_iter1 = last_iter_combiner(last_iters1)
        last_iter2 = last_iter_combiner(last_iters2)
        all_any_checker = ConvergenceChecker(;
            norm_condition = all_condition,
            component_condition = any_condition,
            condition_combiner,
        )
        any_all_checker = ConvergenceChecker(;
            norm_condition = any_condition,
            component_condition = all_condition,
            condition_combiner,
        )
        for checker in (all_any_checker, any_all_checker)
            @test test_func(checker, last_iter1, last_iter2)
        end
    end
end
