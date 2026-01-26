module MultiBroadcastFusionCUDAExt

import CUDA, Adapt
import MultiBroadcastFusion as MBF
import MultiBroadcastFusion: fused_copyto!

MBF.device(x::CUDA.CuArray) = MBF.MBF_CUDA()

include("parameter_memory.jl")

"""
    partition_kernels(fmb;
        fused_broadcast_constructor = MBF.FusedMultiBroadcast,
        args_func::Function = 
    )

Splits fused broadcast kernels into a vector
of kernels, based on parameter memory limitations.

We first attempt to fuse
    1:N, 1:N-1, 1:N-2, ... until we fuse 1:N-k
Next, we attempt to fuse
    N-k+1:N, N-k+1:N-1, N-k+1:N-2, ...

And so forth.
"""
function partition_kernels(
    fmb,
    fused_broadcast_constructor = MBF.FusedMultiBroadcast,
    args_func::Function = fused_multibroadcast_args,
)
    plim = get_param_lim()
    usage = param_usage_args(args_func(fmb))
    n_bins = 1
    fmbs = (fmb,)
    usage ≤ plim && return fmbs
    fmbs_split = []
    N = length(fmb.pairs)
    i_start = 1
    i_stop = N
    while i_stop ≠ i_start
        ith_pairs = fmb.pairs[i_start:i_stop]
        ith_fmb = fused_broadcast_constructor(ith_pairs)
        if param_usage_args(args_func(ith_fmb)) ≤ plim # first iteration will likely fail (ambitious)
            push!(fmbs_split, ith_fmb)
            i_stop == N && break
            i_start = i_stop + 1 # N on first iteration
            i_stop = N # reset i_stop
        else
            i_stop = i_stop - 1
        end
    end
    return fmbs_split
end

function fused_copyto!(fmb::MBF.FusedMultiBroadcast, ::MBF.MBF_CUDA)
    destinations = map(p -> p.first, fmb.pairs)
    fmbs = partition_kernels(fmb)
    for fmb in fmbs
        (; pairs) = fmb
        dest = first(pairs).first
        dests = map(p -> p.first, pairs)
        all(a -> axes(a) == axes(dest), dests) || error(
            "Cannot fuse broadcast expressions with unequal broadcast axes",
        )
        nitems = length(parent(dest))
        CI = CartesianIndices(axes(dest))
        kernel =
            CUDA.@cuda always_inline = true launch = false fused_copyto_kernel!(
                fmb,
                CI,
            )
        config = CUDA.launch_configuration(kernel.fun)
        threads = min(nitems, config.threads)
        blocks = cld(nitems, threads)
        kernel(fmb, CI; threads, blocks)
    end
    return destinations
end
import Base.Broadcast
function fused_copyto_kernel!(fmb::MBF.FusedMultiBroadcast, CI)
    @inbounds begin
        (; pairs) = fmb
        dest = first(pairs).first
        nitems = length(dest)
        idx =
            CUDA.threadIdx().x +
            (CUDA.blockIdx().x - Int32(1)) * CUDA.blockDim().x
        if 1 ≤ idx ≤ nitems
            MBF.rcopyto_at!(pairs, CI[idx])
        end
    end
    return nothing
end

adapt_f(to, f::F) where {F} = Adapt.adapt(to, f)
adapt_f(to, ::Type{F}) where {F} = (x...) -> F(x...)

adapt_src(to, src::AbstractArray) = Adapt.adapt(to, src)

function adapt_src(to, bc::Base.Broadcast.Broadcasted)
    Base.Broadcast.Broadcasted(
        bc.style,
        adapt_f(to, bc.f),
        Adapt.adapt(to, bc.args),
        Adapt.adapt(to, bc.axes),
    )
end

function Adapt.adapt_structure(
    to::CUDA.KernelAdaptor,
    fmbc::MBF.FusedMultiBroadcast,
)
    MBF.FusedMultiBroadcast(map(fmbc.pairs) do pair
        dest = pair.first
        src = pair.second
        Pair(Adapt.adapt(to, dest), adapt_src(to, src))
    end)
end

end
