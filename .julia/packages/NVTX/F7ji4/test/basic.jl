using NVTX, Test

using Colors

@assert NVTX.isactive()

domain = NVTX.Domain("Custom domain")
NVTX.name_category(domain, 1, "boo")
NVTX.name_category(domain, 2, "blah")


NVTX.mark(domain; message="mark 1", category=1)

outer_range = NVTX.range_start(domain; message="outer range", color=colorant"red")

Threads.@threads for i = 1:5
    NVTX.range_push(domain; message="inner range", category=2, payload=i)
    sleep(0.1)
    NVTX.range_pop(domain)
end

NVTX.range_end(outer_range)

GC.gc(false)
GC.gc(true)

module TestMod
using NVTX

const cat = NVTX.@category 34 "my category"

function dostuff(x)
    NVTX.@mark "a mark" category=cat
    NVTX.@mark "mark $x" payload=x
    NVTX.@range x += 1
    NVTX.@range "sleeping" begin
        y = true
        sleep(0.3)
    end
    return y
end

NVTX.@annotate function dostuff(x::Float64)
    sleep(x)
    return true
end

end

@test NVTX.Domain(TestMod) === NVTX.Domain(TestMod)

@test TestMod.dostuff(1)
@test TestMod.dostuff(2)
@test TestMod.dostuff(0.3)
