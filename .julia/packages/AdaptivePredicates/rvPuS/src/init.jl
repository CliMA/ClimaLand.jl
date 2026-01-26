struct InitConstants{T}
    epsilon::T
    splitter::T
    resulterrbound::T
    ccwerrboundA::T
    ccwerrboundB::T
    ccwerrboundC::T
    o3derrboundA::T
    o3derrboundB::T
    o3derrboundC::T
    iccerrboundA::T
    iccerrboundB::T
    iccerrboundC::T
    isperrboundA::T
    isperrboundB::T
    isperrboundC::T
end

function InitConstants{T}() where {T}
    every_other = true
    half = one(T) / 2
    epsilon = one(T)
    splitter = one(T)
    check = one(T)
    cond = true
    while cond
        lastcheck = check
        epsilon *= half
        if every_other
            splitter *= 2one(T)
        end
        every_other = !every_other
        check = one(T) + epsilon
        cond = !isone(check) && check ≠ lastcheck
    end
    @assert epsilon == eps(T) / 2
    splitter += one(T)
    resulterrbound = (3 + 8epsilon) * epsilon
    ccwerrboundA = (3 + 16epsilon) * epsilon
    ccwerrboundB = (2 + 12 * epsilon) * epsilon
    ccwerrboundC = (9 + 64epsilon) * epsilon^2
    o3derrboundA = (7 + 56epsilon) * epsilon
    o3derrboundB = (3 + 28epsilon) * epsilon
    o3derrboundC = (26 + 288epsilon) * epsilon^2
    iccerrboundA = (10 + 96epsilon) * epsilon
    iccerrboundB = (4 + 48epsilon) * epsilon
    iccerrboundC = (44 + 576epsilon) * epsilon^2
    isperrboundA = (16 + 224epsilon) * epsilon
    isperrboundB = (5 + 72epsilon) * epsilon
    isperrboundC = (71 + 1408epsilon) * epsilon^2
    return InitConstants{T}(
        epsilon,
        splitter,
        resulterrbound,
        ccwerrboundA,
        ccwerrboundB,
        ccwerrboundC,
        o3derrboundA,
        o3derrboundB,
        o3derrboundC,
        iccerrboundA,
        iccerrboundB,
        iccerrboundC,
        isperrboundA,
        isperrboundB,
        isperrboundC
    )
end

const IC64 = InitConstants{Float64}()
const IC32 = InitConstants{Float32}()
for e in fieldnames(InitConstants)
    @eval @inline $e(::Type{Float64}) = IC64.$e 
    @eval @inline $e(::Type{Float32}) = IC32.$e 
end