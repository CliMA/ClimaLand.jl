using CpuId

num_machines() = static(1)
num_sockets() = static(1)

_get_num_cores()::Int = clamp(CpuId.cpucores(), 1, (get_cpu_threads())::Int)

const nc = @load_preference("nc", _get_num_cores())
const syst = @load_preference("syst", get_cpu_threads())

num_l1cache() = static(nc)
num_cores() = static(nc)
sys_threads() = static(syst)
num_l2cache() = num_l1cache()
num_l3cache() = static(1)
num_l4cache() = static(0)
cache_inclusive(_) = False()

# TODO: implement
cache_associativity(_) = static(0)

cache_type(::Union{Val{1},StaticInt{1}}) = Val{:Data}()
cache_type(_) = Val{:Unified}()
# cache_type(::Union{Val{2},StaticInt{2}}) = Val{:Unified}()
# cache_type(::Union{Val{3},StaticInt{3}}) = Val{:Unified}()
let lnsize = static(CpuId.cachelinesize())
  global cache_linesize(_) = lnsize
end
cache_size(_) = StaticInt{0}()

const cs = @load_preference("cs", 
  let cs = CpuId.cachesize()
    ntuple(i -> i == 3 ? cs[3] ÷ _get_num_cores() : cs[i], length(cs))
  end)
const ci = @load_preference("ci", CpuId.cacheinclusive())

for (i, csi) in enumerate(cs)
  @eval cache_size(::Union{Val{$i},StaticInt{$i}}) = $(static(csi))
end
for (i, cii) in enumerate(ci)
  @eval cache_inclusive(::Union{Val{$i},StaticInt{$i}}) = $(static(cii != 0))
end
