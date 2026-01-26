# Utilities for help in accessing/creating the relevant directories.

module User

using ..BaseDirs
using ..Internals

@static if VERSION >= v"1.11"
    eval(Expr(:public, :bin, :data, :config, :state, :cache, :runtime,
              :fonts, :applications, :desktop, :downloads, :documents,
              :music, :pictures, :videos, :templates, :public))
end

function bin(; create::Bool=false)
    path = Internals.resolvedirpath(BaseDirs.BIN_HOME, (); create)
    create && Internals.ensureexecutable(path)
    path
end
function bin(name; create::Bool=false)
    path = Internals.resolvedirpath(BaseDirs.BIN_HOME, (name,); create)
    create && Internals.ensureexecutable(path)
    path
end
function bin(project::BaseDirs.Project; create::Bool=false)
    path = Internals.resolvedirpath(BaseDirs.BIN_HOME, (project.name,); create)
    create && Internals.ensureexecutable(path)
    path
end

@defaccessor data DATA_HOME
@defaccessor config CONFIG_HOME
@defaccessor state STATE_HOME
@defaccessor cache CACHE_HOME
@defaccessor runtime RUNTIME_DIR

fonts(pathcomponents...; kwargs...) =
    Internals.resolvedirpaths(filter(p -> startswith(p, homedir()), BaseDirs.FONTS_DIRS),
                              pathcomponents; kwargs...)
applications(pathcomponents...; kwargs...) =
    Internals.resolvedirpaths(filter(p -> startswith(p, homedir()), BaseDirs.APPLICATIONS_DIRS),
                              pathcomponents; kwargs...)

desktop(pathcomponents...)   = joinpath(BaseDirs.DESKTOP_DIR,     pathcomponents...)
downloads(pathcomponents...) = joinpath(BaseDirs.DOWNLOAD_DIR,    pathcomponents...)
documents(pathcomponents...) = joinpath(BaseDirs.DOCUMENTS_DIR,   pathcomponents...)
music(pathcomponents...)     = joinpath(BaseDirs.MUSIC_DIR,       pathcomponents...)
pictures(pathcomponents...)  = joinpath(BaseDirs.PICTURES_DIR,    pathcomponents...)
videos(pathcomponents...)    = joinpath(BaseDirs.VIDEOS_DIR,      pathcomponents...)
templates(pathcomponents...) = joinpath(BaseDirs.TEMPLATES_DIR,   pathcomponents...)
public(pathcomponents...)    = joinpath(BaseDirs.PUBLICSHARE_DIR, pathcomponents...)

end

# ---------

module System

using ..BaseDirs
using ..Internals

@static if VERSION >= v"1.11"
    eval(Expr(:public, :data, :config, :fonts, :applications))
end

@defaccessor data DATA_DIRS
@defaccessor config CONFIG_DIRS

fonts(pathcomponents...; kwargs...) =
    Internals.resolvedirpaths(filter(p -> !startswith(p, homedir()), BaseDirs.FONTS_DIRS),
                              pathcomponents; kwargs...)
applications(pathcomponents...; kwargs...) =
    Internals.resolvedirpaths(filter(p -> !startswith(p, homedir()), BaseDirs.APPLICATIONS_DIRS),
                              pathcomponents; kwargs...)

end

# ---------

@defaccessor data vcat(DATA_HOME, DATA_DIRS)
@defaccessor config vcat(CONFIG_HOME, CONFIG_DIRS)
@defaccessor fonts FONTS_DIRS
@defaccessor applications APPLICATIONS_DIRS

# User accessors without system equivalents, for convenience.
const cache = User.cache
const runtime = User.runtime
const state = User.state
