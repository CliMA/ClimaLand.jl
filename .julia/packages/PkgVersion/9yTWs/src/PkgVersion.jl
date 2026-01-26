module PkgVersion
using Pkg

function project_data(m::Module, name, T, default)
    function _pkgdir(m::Module)
        @static if Base.VERSION < v"1.4"
            pf = Base.pathof(Base.moduleroot(m))
            pf === nothing ? nothing : abspath(pf, "..", "..")
        else
            pkgdir(m)
        end
    end
    pf = _pkgdir(m)
    pf === nothing && return T(default)
    pf = Pkg.Types.projectfile_path(pf)
    project_data(pf, name, T, default)
end

function project_variable(pr::Dict, name, T, default)
    res = get(() -> T(default), pr, string(name))
    res = res isa AbstractVector ? res[1] : res
    T isa Type ? T(res) : res
end
function project_variable(pr, name, T, default)
    sname = Symbol(name)
    res = isdefined(pr, sname) ? getfield(pr, sname) : nothing
    res === nothing ? project_variable(pr.other, name, T, default) : res
end

function project_data(pf::AbstractString, name, T::Union{Type,Function}, default)
    pf === nothing && return T(default)
    pr = Pkg.Operations.read_project(pf)
    project_variable(pr, name, T, default)
end

macro Version(default=0)
    project_data(__module__, :version, VersionNumber, default)
end
macro Uuid(default=0)
    project_data(__module__, :uuid, Base.UUID, default)
end
macro Author(default="unknown")
    project_data(__module__, "authors", string, default)
end

Version(m::Module, default=0) = project_data(m, :version, VersionNumber, default)
Uuid(m::Module, default=0) = project_data(m, :uuid, Base.UUID, default)
Author(m::Module, default="unknown") = project_data(m, "authors", string, default)

const VERSION = PkgVersion.@Version

end # module
