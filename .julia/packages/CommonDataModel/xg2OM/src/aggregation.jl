
aggregator(::Type{typeof(var)}) = VarianceWelfordAggregation
aggregator(::Type{typeof(maximum)}) = MaxAggregation
aggregator(::Type{typeof(minimum)}) = MinAggregation

struct VarianceWelfordAggregation{Ti,T}
    count::Ti
    mean::T
    M2::T
end

function VarianceWelfordAggregation(T)
    VarianceWelfordAggregation{Int,T}(0,zero(T),zero(T))
end

# Welford's online algorithm
@inline function update(ag::VarianceWelfordAggregation{Ti,T}, new_value) where {Ti,T}
    count = ag.count
    mean = ag.mean
    M2 = ag.M2

    count += 1
    delta = new_value - mean
    mean += delta / count
    delta2 = new_value - mean
    M2 += delta * delta2
    return VarianceWelfordAggregation{Ti,T}(count, mean, M2)
end

Base.count(ag::VarianceWelfordAggregation) = ag.count

function Statistics.mean(ag::VarianceWelfordAggregation{Ti,T}) where {Ti,T}
    if ag.count > 0
        return ag.mean
    else
        # NaN of the right kind
        return zero(T) / 0
    end
end

function Statistics.var(ag::VarianceWelfordAggregation{Ti,T}) where {Ti,T}
    if ag.count > 1
        sample_variance = ag.M2 / (ag.count - 1)
    else
        # ag.count is 0 or 1
        # NaN of the right kind
        sample_variance = zero(T) / 0
    end
    return sample_variance
end

result(ag::VarianceWelfordAggregation) = var(ag)


for (funAggregation,fun) in ((:MaxAggregation,max),(:MinAggregation,min))
    @eval begin
        struct $funAggregation{T}
            result::T
            init::Bool
        end

        function $funAggregation(T)
            $funAggregation{T}(zero(T),false)
        end

        @inline function update(ag::$funAggregation{T}, new_value) where T
            if ag.init
                return $funAggregation{T}(max(ag.result,new_value),true)
            else
                return $funAggregation{T}(new_value,true)
            end
        end

        function result(ag::$funAggregation)
            if ag.init
                return ag.result
            else
                error("reducing over an empty collection is not allowed")
            end
        end
    end
end
