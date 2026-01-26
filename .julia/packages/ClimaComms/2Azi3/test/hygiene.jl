import ClimaComms as CC
CC.@import_required_backends
function test_macro_hyhiene(dev)
    AT = CC.array_type(dev)
    n = 3 # tests that we can reach variables defined in scope

    array = AT(rand(10))
    CC.@threaded dev for i in 1:n
        array[i] = sin(array[i])
    end

    CC.time(dev) do
        for i in 1:n
            sin.(AT(rand(1000, 1000)))
        end
    end

    CC.@time dev for i in 1:n
        sin.(AT(rand(1000, 1000)))
    end

    CC.elapsed(dev) do
        for i in 1:n
            sin.(AT(rand(10)))
        end
    end

    CC.@elapsed dev for i in 1:n
        sin.(AT(rand(10)))
    end

    CC.@assert dev true
    CC.@assert dev true "some message: $(true)"

    CC.sync(dev) do
        for i in 1:n
            sin.(AT(rand(10)))
        end
    end

    CC.@sync dev for i in 1:n
        sin.(AT(rand(10)))
    end

    CC.cuda_sync(dev) do
        for i in 1:n
            sin.(AT(rand(10)))
        end
    end

    CC.@cuda_sync dev for i in 1:n
        sin.(AT(rand(10)))
    end
end
dev = CC.device()

test_macro_hyhiene(dev)
