using Hwloc

mutable struct Topology
  topology::Union{Nothing,Hwloc.Object}
end
Topology() = Topology(nothing)
function safe_topology_load!()
  try
    TOPOLOGY.topology = Hwloc.gettopology()
  catch e
    @warn e
    @warn """
          Using Hwloc failed. Please file an issue with the above warning at: https://github.com/JuliaParallel/Hwloc.jl
          Proceeding with generic topology assumptions. This may result in reduced performance.
      """
  end
end

const TOPOLOGY = Topology()
safe_topology_load!()

function count_attr(topology::Hwloc.Object, attr::Symbol)
  count::Int = 0
  for t ∈ topology
    count += t.type_ === attr
  end
  if ((Sys.ARCH === :aarch64) && Sys.isapple()) && ((attr === :Core) | (attr === :PU))
    # count big cores only on Apple, as one process can only use big or small cores at a time.
    count >>>= one(count) # FIXME: how to actually identify big cores???
  end
  count
end
function count_attr(attr::Symbol)
  topology = TOPOLOGY.topology
  topology === nothing && return nothing
  count_attr(topology, attr)
end

@static if Sys.isapple() && Sys.ARCH === :aarch64 # detect M1
  redefine_attr_count() = nothing
  num_l1cache() = static(4)
  num_l2cache() = static(1)
  num_l3cache() = static(0)
  num_l4cache() = static(0)
  num_machines() = static(1)
  num_sockets() = static(1)
  num_cores() = static(4)
  sys_threads() = static(4)
else # not M1
  @noinline function define_attr_count(fname::Symbol, v)
    if v === nothing
      @eval $fname() = nothing
      return
    elseif v isa Integer
      @eval $fname() = StaticInt{$(convert(Int, v))}()
    elseif v isa Bool
      if v
        @eval $fname() = True()
      else
        @eval $fname() = False()
      end
    else
      @eval $fname() = $v
    end
    nothing
  end
  # define_attr(fname::Symbol, attr::Symbol) = define_attr_count(fname, )

  for (f, attr) ∈ [
    (:num_l1cache, :L1Cache),
    (:num_l2cache, :L2Cache),
    (:num_l3cache, :L3Cache),
    (:num_l4cache, :L4Cache),
    (:num_machines, :Machine),
    (:num_sockets, :Package),
    (:num_cores, :Core),
    (:sys_threads, :PU),
  ]
    define_attr_count(f, count_attr(attr))
  end

  function redefine_attr_count()
    sys_thread = Int(sys_threads())::Int
    iter = [
      (num_l1cache(), :num_l1cache, :L1Cache),
      (num_l2cache(), :num_l2cache, :L2Cache),
      (num_l3cache(), :num_l3cache, :L3Cache),
      (num_l4cache(), :num_l4cache, :L4Cache),
      (num_machines(), :num_machines, :Machine),
      (num_sockets(), :num_sockets, :Package),
      (num_cores(), :num_cores, :Core),
      (sys_thread, :sys_threads, :PU),
    ]
    for (v, f, attr) ∈ iter
      ref = count_attr(attr)
      if ref ≠ v
        @debug "Redefining attr count $f = $ref."
        define_attr_count(f, ref)
      end
    end
    nothing
  end
end # not M1

function dynamic_cache_inclusivity()::NTuple{4,Bool}
  @static if !((Sys.ARCH === :x86_64) || (Sys.ARCH === :i686))
    return (false, false, false, false)
  end
  function get_cache_edx(subleaf)
    # source: https://github.com/m-j-w/CpuId.jl/blob/401b638cb5a020557bce7daaf130963fb9c915f0/src/CpuInstructions.jl#L38
    # credit Markus J. Weber, copyright: https://github.com/m-j-w/CpuId.jl/blob/master/LICENSE.md
    Base.llvmcall(
      """
              ; leaf = %0, subleaf = %1, %2 is some label
              ; call 'cpuid' with arguments loaded into registers EAX = leaf, ECX = subleaf
              %2 = tail call { i32, i32, i32, i32 } asm sideeffect "cpuid",
                  "={ax},={bx},={cx},={dx},{ax},{cx},~{dirflag},~{fpsr},~{flags}"
                  (i32 4, i32 %0) #2
              ; retrieve the result values and return eax and edx contents
              %3 = extractvalue { i32, i32, i32, i32 } %2, 0
              %4 = extractvalue { i32, i32, i32, i32 } %2, 3
              %5  = insertvalue [2 x i32] undef, i32 %3, 0
              %6  = insertvalue [2 x i32]   %5 , i32 %4, 1
              ; return the value
              ret [2 x i32] %6
              """
      # llvmcall requires actual types, rather than the usual (...) tuple
      ,
      Tuple{UInt32,UInt32},
      Tuple{UInt32},
      subleaf % UInt32,
    )
  end
  # eax0, edx1 = get_cache_edx(0x00000000)
  t = (false, false, false, false)
  i = zero(UInt32)
  j = 0
  while (j < 4)
    eax, edx = get_cache_edx(i)
    i += one(UInt32)
    iszero(eax & 0x1f) && break
    iszero(eax & 0x01) && continue
    ci = ((edx & 0x00000002) != 0x00000000) & (eax & 0x1f != 0x00000000)
    t = Base.setindex(t, ci, (j += 1))
  end
  t
end

nothing_cache_summary() =
  (size = 0, linesize = 64, associativity = nothing, type = nothing, inclusive = nothing)
function dynamic_cache_summary(N)
  topology = TOPOLOGY.topology
  cache_name = (:L1Cache, :L2Cache, :L3Cache, :L4Cache)[N]
  if (topology === nothing) || (iszero(count_attr(cache_name)))
    return nothing_cache_summary()
  end
  c = first(t for t in topology if t.type_ == cache_name && t.attr.depth == N).attr
  (
    size = c.size,
    linesize = c.linesize,
    associativity = c.associativity,
    type = c.type_,
    inclusive = dynamic_cache_inclusivity()[N],
  )
end
cache_size(_) = StaticInt{0}()
cache_associativity(_) = nothing
cache_type(_) = nothing
cache_inclusive(_) = nothing

unwrap(::Val{S}) where {S} = S
@static if Sys.isapple() && Sys.ARCH === :aarch64
  cache_linesize(_) = StaticInt{128}() # assume...
  redefine_cache(_) = nothing
  cache_size(::Union{Val{1},StaticInt{1}}) = StaticInt{131072}()
  cache_linesize(::Union{Val{1},StaticInt{1}}) = StaticInt{128}()
  cache_associativity(::Union{Val{1},StaticInt{1}}) = StaticInt{0}()
  cache_type(::Union{Val{1},StaticInt{1}}) = Val{:Data}()
  cache_inclusive(::Union{Val{1},StaticInt{1}}) = False()

  cache_size(::Union{Val{2},StaticInt{2}}) = StaticInt{12582912}()
  cache_linesize(::Union{Val{2},StaticInt{2}}) = StaticInt{128}()
  cache_associativity(::Union{Val{2},StaticInt{2}}) = StaticInt{0}()
  cache_type(::Union{Val{2},StaticInt{2}}) = Val{:Unified}()
  cache_inclusive(::Union{Val{2},StaticInt{2}}) = False()

else # not M1
  cache_linesize(_) = StaticInt{64}() # assume...
  function define_cache(N, c = dynamic_cache_summary(N))
    c === nothing_cache_summary() || _define_cache(N, c)
  end
  @noinline function _define_cache(N, c)
    @eval begin
      cache_size(::Union{Val{$N},StaticInt{$N}}) = StaticInt{$(c.size)}()
      cache_linesize(::Union{Val{$N},StaticInt{$N}}) = StaticInt{$(c.linesize)}()
      cache_associativity(::Union{Val{$N},StaticInt{$N}}) = StaticInt{$(c.associativity)}()
      cache_type(::Union{Val{$N},StaticInt{$N}}) =
        Val{$(c.type === nothing ? nothing : QuoteNode(c.type))}()
      cache_inclusive(::Union{Val{$N},StaticInt{$N}}) =
        $(c.inclusive isa Bool && c.inclusive ? :True : :False)()
    end
    nothing
  end
  function redefine_cache(N)
    s = cache_size(StaticInt(N))
    l = cache_linesize(StaticInt(N))
    a = cache_associativity(StaticInt(N))
    t = cache_type(StaticInt(N))
    i = cache_inclusive(StaticInt(N))
    c = (
      size = Int(s)::Int,
      linesize = l === nothing ? nothing : Int(l)::Int,
      associativity = a === nothing ? nothing : Int(a)::Int,
      type = t === nothing ? nothing : unwrap(t)::Symbol,
      inclusive = i === nothing ? nothing : Bool(i)::Bool,
    )
    correct = dynamic_cache_summary(N)
    if c !== correct
      @debug "Redefining cache $N."
      _define_cache(N, correct)
    end
    nothing
  end
  foreach(define_cache, 1:4)
end

function __init__()
  Sys.isapple() && Sys.ARCH === :aarch64 && return # detect M1
  ccall(:jl_generating_output, Cint, ()) == 1 && return
  safe_topology_load!()
  if count_attr(:Core) ≢ BASELINE_CORES
    redefine_attr_count()
    foreach(redefine_cache, 1:4)
  end
  return nothing
end
