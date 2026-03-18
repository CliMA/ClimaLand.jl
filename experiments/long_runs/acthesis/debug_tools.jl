@inline join_names(headname::String, propertyname) =
    headname * "." * string(propertyname)

@inline get_tendency_object(::Val{ClimaTimeSteppers.ARS111()}, simulation) =
    simulation._integrator.cache.T_exp[1]

function check_nans(
    x::ClimaCore.Fields.FieldVector,
    mask,
    pname::String = "X";
    excl = (),
    verbose::Bool = false,
)
    return any(
        check_nans(
            getproperty(x, p),
            mask,
            join_names(pname, p),
            excl = excl,
            verbose = verbose,
        ) for p in propertynames(x) if !(p in excl)
    )
end

function check_nans(
    x::NamedTuple,
    mask,
    pname::String = "X";
    excl = (),
    verbose::Bool = false,
)
    return any(
        check_nans(
            x[p],
            mask,
            join_names(pname, p),
            excl = excl,
            verbose = verbose,
        ) for p in propertynames(x) if !(p in excl)
    )
end

function check_nans(
    x::ClimaCore.Fields.Field,
    mask,
    pname::String = "X";
    excl = (),
    verbose::Bool = false,
)
    return check_nans(
        x,
        propertynames(x),
        mask,
        pname,
        excl = excl,
        verbose = verbose,
    )
end

function check_nans(
    x::ClimaCore.Fields.Field,
    pnames::Tuple{Vararg{Any}},
    mask,
    headname::String;
    excl = (),
    verbose::Bool = false,
)
    return any(
        check_nans(
            getproperty(x, p),
            mask,
            join_names(headname, p),
            excl = excl,
            verbose = verbose,
        ) for p in pnames if !(p in excl)
    )
end

function check_nans(
    x::ClimaCore.Fields.Field,
    pnames::Tuple{},
    mask,
    pname::String;
    excl = (),
    verbose::Bool = false,
)
    return check_nans(x, axes(x), mask, pname; verbose = verbose)
end

function check_nans(
    x::ClimaCore.Fields.Field,
    space::Union{
        ClimaCore.Spaces.AbstractSpectralElementSpace,
        ClimaCore.Spaces.AbstractPointSpace,
    },
    mask,
    pname::String;
    verbose::Bool = false,
)
    return count_nans(x, mask, pname, Val(verbose))
end

function check_nans(
    x::ClimaCore.Fields.Field,
    space::Union{
        ClimaCore.Spaces.ExtrudedFiniteDifferenceSpace,
        ClimaCore.Spaces.FiniteDifferenceSpace,
    },
    mask,
    pname::String;
    verbose::Bool = false,
)
    return count_nans(
        ClimaLand.Domains.top_center_to_surface(x),
        mask,
        pname,
        Val(verbose),
    )
end

function count_nans(
    x::ClimaCore.Fields.Field,
    mask,
    pname::String,
    verbose::Val{false},
)
    isnothing(mask) && return count(isnan, parent(x)) > 0
    return mapreduce((s, m) -> m != 0 && isnan(s), |, parent(x), parent(mask))
end

function count_nans(
    x::ClimaCore.Fields.Field,
    mask,
    pname::String,
    verbose::Val{true},
)
    n =
        isnothing(mask) ? count(isnan, parent(x)) :
        mapreduce(
            (s, m) -> m != 0 && isnan(s),
            Base.add_sum,
            parent(x),
            parent(mask),
        )
    print(pname, ": ", n, "\n")
    return n >= 1
end

function nans_exist(simulation, mask; excl = (), verbose::Bool = false)
    nans_in_Y = check_nans(
        simulation._integrator.u,
        mask,
        "Y",
        excl = excl,
        verbose = verbose,
    )
    nans_in_p = check_nans(
        simulation._integrator.p,
        mask,
        "p",
        excl = excl,
        verbose = verbose,
    )
    nans_in_dY = check_nans(
        get_tendency_object(Val(simulation.timestepper.name), simulation),
        mask,
        "dY",
        excl = excl,
        verbose = verbose,
    )
    return nans_in_Y || nans_in_p || nans_in_dY
end

function write_present_state(simulation, old_Y, old_p, outdir)
    curr_dY = get_tendency_object(Val(simulation.timestepper.name), simulation)
    curr_Y = simulation._integrator.u
    curr_p = simulation._integrator.p
    JLD2.jldsave(
        joinpath(outdir, "integrator_state.jld2"),
        curr_t = cpu(simulation._integrator.t),
        Y = my_to_cpu(curr_Y),
        p = my_to_cpu(curr_p),
        dY = my_to_cpu(curr_dY),
        old_Y = my_to_cpu(old_Y),
        old_p = my_to_cpu(old_p),
    )
end

function load_present_state(filepath)
    data = JLD2.load(filepath)
    return NamedTuple(Symbol(k) => v for (k, v) in data)
end

function any_negatives(x::ClimaCore.Fields.Field, mask)
    isnothing(mask) && return count(<(0), parent(x)) > 0
    return mapreduce((s, m) -> m != 0 && s < 0, |, parent(x), parent(mask))
end

########################################

function set_fields!(
    dest::ClimaCore.Fields.FieldVector,
    x::ClimaCore.Fields.FieldVector,
)
    for p in propertynames(x)
        set_fields!(getproperty(dest, p), getproperty(x, p))
    end
end

function set_fields!(dest::NamedTuple, x::NamedTuple)
    for p in propertynames(x)
        set_fields!(dest[p], x[p])
    end
end

function set_fields!(dest::ClimaCore.Fields.Field, x::ClimaCore.Fields.Field)
    pnames = propertynames(x)
    if pnames isa Tuple{}
        @. dest = x
    else
        for p in pnames
            set_fields!(getproperty(dest, p), getproperty(x, p))
        end
    end
end

my_to_cpu(x) = ClimaCore.to_cpu(x)
my_to_cpu(x::NamedTuple) = map(my_to_cpu, x)
#my_to_gpu(x) = ClimaCore.to_device(device, x)
#my_to_gpu(x::NamedTuple) = map(my_to_gpu, x)

using NCDatasets
function get_bad_data(root_dir; use_idx = 0)
    zdata = NCDataset(joinpath(root_dir, "snd_1d_average.nc"))
    swedata = NCDataset(joinpath(root_dir, "swe_1d_average.nc"))
    swadata = NCDataset(joinpath(root_dir, "swa_1d_average.nc"))
    scfdata = NCDataset(joinpath(root_dir, "snowc_1d_average.nc"))
    galbdata = NCDataset(joinpath(root_dir, "galb_1d_average.nc"))
    salbdata = NCDataset(joinpath(root_dir, "salb_1d_average.nc"))
    snalbdata = NCDataset(joinpath(root_dir, "snalb_1d_average.nc"))
    swe = swedata[:swe][:, :, :]
    tidx = use_idx != 0 ? use_idx : length(swe[:, 1, 1])
    set2 = findall(isnan, swedata[:swe][tidx, :, :])
    set1 = findall(isnan, swedata[:swe][1, :, :])
    idxs = collect(setdiff(set2, set1))
    bad_tij = [(findfirst(isnan, swe[:, x[1], x[2]]), x[1], x[2]) for x in idxs]
    bad_data = [
        (
            t_idx = t,
            lon_idx = i,
            lat_idx = j,
            t = zdata[:date][t],
            lon = zdata[:lon][i],
            lat = zdata[:lat][j],
            z = zdata[:snd][t - 1, i, j],
            swe = swedata[:swe][t - 1, i, j],
            scf = scfdata[:snowc][t - 1, i, j],
            cell_alb = swadata[:swa][t - 1, i, j],
            snow_alb = snalbdata[:snalb][t - 1, i, j],
            soil_alb = salbdata[:salb][t - 1, i, j],
            surf_alb = galbdata[:galb][t - 1, i, j],
        ) for (t, i, j) in bad_tij
    ]
    return bad_data
end

function get_col_timeseries(root_dir, col_ij)
    files_in_dir = readdir(root_dir)
    data = Dict()
    var_data = []
    for file in files_in_dir
        print(file, "\n")
        if file[(end - 2):end] == ".nc"
            var_data = NCDataset(joinpath(root_dir, file))
            var_name = Symbol(split(file, "_")[1])
            if var_name != :tsoil
                data[var_name] = var_data[var_name][:, col_ij[1], col_ij[2]]
            else
                data[var_name] = var_data[var_name][:, col_ij[1], col_ij[2], :]
            end
        end
    end
    data[:t] = var_data[:date][:]
    return NamedTuple(data)
end
