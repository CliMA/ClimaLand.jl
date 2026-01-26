using ClimaComms

ClimaComms.@import_required_backends

context = ClimaComms.context()
pid, nprocs = ClimaComms.init(context)
device = ClimaComms.device(context)
AT = ClimaComms.array_type(device)

x_max = 100
y_max = device isa ClimaComms.CUDADevice ? 1000000 : 10000
a = AT(rand(x_max, y_max))
∂ⁿa∂xⁿ = similar(a)

Base.@propagate_inbounds function nth_deriv_along_axis1(array, n, i, indices...)
    @assert n >= 0
    n == 0 && return array[i, indices...]
    prev_i = i == 1 ? size(array, 1) : i - 1
    next_i = i == size(array, 1) ? 1 : i + 1
    value_at_prev_i = nth_deriv_along_axis1(array, n - 1, prev_i, indices...)
    value_at_next_i = nth_deriv_along_axis1(array, n - 1, next_i, indices...)
    value_at_i = nth_deriv_along_axis1(array, n - 1, i, indices...)
    return (value_at_next_i - 2 * value_at_i + value_at_prev_i) / 2
end

function print_time_and_bandwidth(device, reads_and_writes, time)
    @info "    Time = $(round(time; sigdigits = 3)) s"
    if device isa ClimaComms.CUDADevice
        cuda_device = CUDA.device()
        memory_bits_attr = CUDA.DEVICE_ATTRIBUTE_GLOBAL_MEMORY_BUS_WIDTH
        memory_kilohertz_attr = CUDA.DEVICE_ATTRIBUTE_MEMORY_CLOCK_RATE
        peak_bandwidth_gbps =
            (CUDA.attribute(cuda_device, memory_bits_attr) / (8 * 1024^3)) *
            (CUDA.attribute(cuda_device, memory_kilohertz_attr) * 1000)
        actual_bandwidth_gbps = (reads_and_writes / 1024^3) / time
        bandwidth_percent = 100 * actual_bandwidth_gbps / peak_bandwidth_gbps
        @info "    Bandwidth = $(round(bandwidth_percent; sigdigits = 3))% of \
                   $(round(peak_bandwidth_gbps; sigdigits = 3)) GB/s"
    end
end

@info "Benchmarking n-th derivative along first axis of a $x_max×$y_max matrix"

identity_copy_reference!(∂ⁿa∂xⁿ, a, device) =
    ClimaComms.@cuda_sync device ∂ⁿa∂xⁿ .= a
ident_first_time = @elapsed identity_copy_reference!(∂ⁿa∂xⁿ, a, device)
ident_time = @elapsed identity_copy_reference!(∂ⁿa∂xⁿ, a, device)
@info "reference identity copy, n = 0:"
@info "    Latency = $(round(ident_first_time - ident_time; sigdigits = 3)) s"
print_time_and_bandwidth(device, 2 * sizeof(a), ident_time)

identity_copy!(∂ⁿa∂xⁿ, a, device) = ClimaComms.@cuda_sync device begin
    ClimaComms.@threaded device for y in axes(a, 2), x in axes(a, 1)
        @inbounds ∂ⁿa∂xⁿ[x, y] = a[x, y]
    end
end
ident_first_time = @elapsed identity_copy_reference!(∂ⁿa∂xⁿ, a, device)
ident_time = @elapsed identity_copy_reference!(∂ⁿa∂xⁿ, a, device)
@info "@threaded identity copy, n = 0:"
@info "    Latency = $(round(ident_first_time - ident_time; sigdigits = 3)) s"
print_time_and_bandwidth(device, 2 * sizeof(a), ident_time)

xy_pairs = AT(CartesianIndices(a))
nth_deriv_by_col_reference!(∂ⁿa∂xⁿ, a, device, n, xy_pairs) =
    ClimaComms.@cuda_sync device begin
        ∂ⁿa∂xⁿ .=
            (pair -> nth_deriv_along_axis1(a, n, pair[1], pair[2])).(xy_pairs)
    end
for n in (0, 2, 6)
    if n == 0
        first_time =
            @elapsed nth_deriv_by_col_reference!(∂ⁿa∂xⁿ, a, device, n, xy_pairs)
    end
    time = @elapsed nth_deriv_by_col_reference!(∂ⁿa∂xⁿ, a, device, n, xy_pairs)
    @info "reference derivative (broadcast over matrix of indices), n = $n:"
    n == 0 && @info "    Latency = $(round(first_time - time; sigdigits = 3)) s"
    print_time_and_bandwidth(device, 2 * sizeof(a), time)
end

nth_deriv_by_col!(∂ⁿa∂xⁿ, a, device, n) = ClimaComms.@cuda_sync device begin
    ClimaComms.@threaded device for y in axes(a, 2), x in axes(a, 1)
        @inbounds ∂ⁿa∂xⁿ[x, y] = nth_deriv_along_axis1(a, n, x, y)
    end
end
for n in (0, 2, 6)
    if n == 0
        first_time = @elapsed nth_deriv_by_col!(∂ⁿa∂xⁿ, a, device, n)
    end
    time = @elapsed nth_deriv_by_col!(∂ⁿa∂xⁿ, a, device, n)
    @info "@threaded derivative (no shared memory), n = $n:"
    n == 0 && @info "    Latency = $(round(first_time - time; sigdigits = 3)) s"
    print_time_and_bandwidth(device, 2 * sizeof(a), time)
end

# The following shared memory kernels will be benchmarked in PR #113.
#=
nth_deriv_by_col_shmem_1!(∂ⁿa∂xⁿ, a, device, n, ::Val{x_max}) where {x_max} =
    ClimaComms.@cuda_sync device begin
        ClimaComms.@threaded device begin
            for y in axes(a, 2), x in @interdependent(axes(a, 1))
                T = eltype(a)
                a_col = ClimaComms.static_shared_memory_array(device, T, x_max)
                @inbounds begin
                    ClimaComms.@sync_interdependent a_col[x] = a[x, y]
                    ClimaComms.@sync_interdependent ∂ⁿa∂xⁿ[x, y] =
                        nth_deriv_along_axis1(a_col, n, x)
                end
            end
        end
    end

nth_deriv_by_col_shmem_2!(∂ⁿa∂xⁿ, a, device, n, ::Val{x_max}) where {x_max} =
    ClimaComms.@cuda_sync device begin
        ClimaComms.@threaded device begin
            for y in axes(a, 2), x in @interdependent(axes(a, 1))
                T = eltype(a)
                a_col = ClimaComms.static_shared_memory_array(device, T, x_max)
                ∂ᵐa∂xᵐ_col = ClimaComms.static_shared_memory_array(device, T, x_max)
                m = n ÷ 2
                @inbounds begin
                    ClimaComms.@sync_interdependent a_col[x] = a[x, y]
                    ClimaComms.@sync_interdependent ∂ᵐa∂xᵐ_col[x] =
                        nth_deriv_along_axis1(a_col, m, x)
                    ClimaComms.@sync_interdependent ∂ⁿa∂xⁿ[x, y] =
                        nth_deriv_along_axis1(∂ᵐa∂xᵐ_col, n - m, x)
                end
            end
        end
    end
=#
