# Helper functions and macros

module Internals

using ..BaseDirs
using Base.Docs

export @defaccessor, @setxdg, @setxdgs

@static if Sys.isunix()
    macro setxdg(envvar::Symbol, default)
        quote global $(esc(envvar)) = if haskey(ENV, $("XDG_$envvar")) && !isempty(ENV[$("XDG_$envvar")])
            String(chopsuffix(ENV[$("XDG_$envvar")], Base.Filesystem.path_separator))
        else expanduser($(esc(default))) end
        end
    end
else
    macro setxdg(envvar::Symbol, default)
        quote global $(esc(envvar)) = if haskey(ENV, $("XDG_$envvar")) && !isempty(ENV[$("XDG_$envvar")])
            String(chopsuffix(ENV[$("XDG_$envvar")], Base.Filesystem.path_separator))
        else $(esc(default)) end
        end
    end
end

macro setxdgs(envvar::Symbol, defaults)
    quote global $(esc(envvar)) = if haskey(ENV, $("XDG_$envvar")) && !isempty(ENV[$("XDG_$envvar")])
        map(split(ENV[$("XDG_$envvar")], ':')) do path
            String(chopsuffix(path, Base.Filesystem.path_separator))
        end
    else $(esc(defaults)) end
    end
end

const NEW_BASEDIR_MODE = 0o700

ensurebasedir(path::String) =
    isdir(path) || mkpath(path, mode=NEW_BASEDIR_MODE)

const DIRECTORY_SUFFIX_FLAG = '/'

"""
    ensurepath(path::String)

Ensure that `path` exists. Should `path` end with `$DIRECTORY_SUFFIX_FLAG` it is
interpreted as a directory. The directory-suffix `$DIRECTORY_SUFFIX_FLAG` is
used on all filesystems for consistency of the API, regardless of the native
path seperator of the host filesystem.
"""
function ensurepath(path::String)
    if !ispath(path)
        if endswith(path, DIRECTORY_SUFFIX_FLAG) || (Sys.iswindows() && endswith(path, "\\"))
            mkpath(path[begin:prevind(path, end)])
        else
            mkpath(dirname(path))
            touch(path)
        end
    end
end

"""
    ensureexecutable(path::String)

Make `path` executable by everybody who can read it. Returns `path`.
"""
function ensureexecutable(path::String)
    if isfile(path)
        basemode = filemode(path)
        uread = basemode & 0o400 > 0
        gread = basemode & 0o040 > 0
        oread = basemode & 0o004 > 0
        xmask = 0o100 * uread + 0o010 * gread + 0o001 * oread
        chmod(path, basemode | xmask)
    end
    path
end

"""
    PROMISED_NO_CONST_ASSIGNMENT

Whether it is promised that no constant assignment will be made to a `BaseDirs`.
"""
PROMISED_NO_CONST_ASSIGNMENT::Bool = false

function warn_if_precompiling()
    if ccall(:jl_generating_output, Cint, ()) != 0 && !PROMISED_NO_CONST_ASSIGNMENT
        @noinline (function ()
                       st = stacktrace(backtrace())
                       any(sf -> sf.func ∈ (:__init__, :_atexit), st) && return
                       naughtyline = findfirst(sf -> sf.file != Symbol(@__FILE__), st)
                       naughtysf = st[something(naughtyline, 1)]
                       @warn """A base directory is being computed during precompilation.
                                This is dangerous, as results depend on the live system configuration.

                                It is recommended that you invoke BaseDirs as required in
                                function bodies, rather than at the top level. Calls are very
                                cheap, and you can always pass the result of a live call around.

                                If you have verified that this call was not assigned to a global constant,
                                you can silence this warning with `BaseDirs.@promise_no_assign`.
                                """ _module="BaseDirs" _file=String(naughtysf.file) _line=naughtysf.line
                   end)()
    end
end

function warn_misplaced_promise()
    if ccall(:jl_generating_output, Cint, ()) == 0
        @noinline (function ()
                       st = stacktrace(backtrace())
                       naughtyline = findfirst(sf -> sf.file != Symbol(@__FILE__), st)
                       naughtysf = st[something(naughtyline, 1)]
                       @warn """A `BaseDirs.@promise_no_assign` declaration has been evalated at runtime.
                                This macro is explicitly intended for use during precompilation, and is likely misplaced.
                                """ _module="BaseDirs" _file=String(naughtysf.file) _line=naughtysf.line
                   end)()
    end
end

function resolvedirpath(basedir::String, pathcomponents::Union{<:Tuple{Vararg{S}}, <:AbstractVector{S}}; create::Bool=false) where {S<:AbstractString}
    create && ensurebasedir(basedir)
    if isempty(pathcomponents)
        basedir
    else
        fullpath = joinpath(basedir, pathcomponents...)
        create && ensurepath(fullpath)
        fullpath
    end
end

function resolvedirpaths(basedirs::Vector{String}, pathcomponents::Union{<:Tuple{Vararg{S}}, <:AbstractVector{S}}; create::Bool=false, existent::Bool=false) where {S<:AbstractString}
    allpaths = [resolvedirpath(bdir, pathcomponents; create) for bdir in basedirs]
    isreachablepath(path::String) = try ispath(path) catch; false end
    if existent
        filter(isreachablepath, allpaths)
    else
        allpaths
    end
end

macro defaccessor(fnname::Symbol, var::Union{Symbol, Expr})
    dirvar = if var isa Symbol
        GlobalRef(BaseDirs, var)
    else esc(var) end
    vecfns = (:vec, :vcat, :filter, :map, :push!, :pushfirst!) # a few that come to mind
    resolver = if (var isa Symbol && getfield(BaseDirs, var) isa Vector{String}) ||
          (var isa Expr && (var.head == :vect ||
                            (var.head == :call && var.args[1] in vecfns)))
        :resolvedirpaths
    else
        :resolvedirpath
    end
    quote
        function $(esc(fnname))(pathcomponents...; kwargs...)
            $warn_if_precompiling()
            $resolver($dirvar, pathcomponents; kwargs...)
        end
        $(esc(fnname))(project::BaseDirs.Project, pathcomponents...; kwargs...) =
            $(esc(fnname))(BaseDirs.projectpath(project, $dirvar), pathcomponents...; kwargs...)
        $(esc(fnname))(mod::Module, pathcomponents...; kwargs...) =
            $(esc(fnname))(BaseDirs.Project(mod), pathcomponents...; kwargs...)
    end
end

function accessordoc(finfo::Union{Symbol, Tuple{String, Symbol}},
                    var::Union{Nothing, Symbol, Vector{Symbol}}=nothing;
                    plural::Bool=if isnothing(var) false
                    elseif var isa Vector true
                    else getfield(BaseDirs, var) isa Vector{String} end,
                    name::String=String(if finfo isa Symbol finfo else last(finfo) end))
    fprefix, fname = if finfo isa Symbol; ("", finfo) else finfo end
    rettype = ifelse(plural, "Vector{String}", "String")
    dirprefix, dirterm = ifelse(plural, ("all", "directories"), ("the", "directory"))
    existentkwarg = ifelse(plural, " - `existent::Bool` (default `false`), filter out paths that do not exist.", "")
    vardoc = if var isa Vector && (dvars = filter(v -> haskey(Docs.meta(BaseDirs), Docs.Binding(BaseDirs, v)), var)) |> !isempty
        "\nThe returned path is determined by the value of $(join(map(v -> "[`BaseDirs.$v`](@ref)", dvars), ", ", ", and ")), which see.\n"
    elseif var isa Symbol && haskey(Docs.meta(BaseDirs), Docs.Binding(BaseDirs, var))
        "\nThe returned path is determined by the value of [`BaseDirs.$var`](@ref), which see.\n"
    else "" end
    kwargs = ifelse(plural, "; create, existent", "; create")
    """
    $fprefix$fname($kwargs) -> $rettype # $dirprefix $name $dirterm
    $fprefix$fname([mod::Module or proj::Project], parts...$kwargs) # a [project-specific] $name path

Locate $dirprefix $name $dirterm, or a path within it.

A project or module can be optionally provided as the first argument, in which
case the returned path is scoped to the project or module.

$vardoc
## Keyword arguments
 - `create::Bool` (default `false`), whether the path should be created if it
   does not exist. Paths ending in `$DIRECTORY_SUFFIX_FLAG` are interpreted as
   directories, and all other paths are considered files. This takes care to
   create the base directories with the appropriate permissions ($(string(NEW_BASEDIR_MODE, base=8))).
$existentkwarg"""
end

end
