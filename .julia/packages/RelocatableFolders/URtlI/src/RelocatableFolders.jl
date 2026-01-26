module RelocatableFolders

import Scratch, SHA

export @path, getroot

function safe_isfile(file)
    return try
        isfile(file)
    catch err
        err isa Base.IOError || rethrow()
        false
    end
end

function safe_ispath(file)
    return try
        ispath(file)
    catch err
        err isa Base.IOError || rethrow()
        false
    end
end

"""
    @path expr [ignore]

Define a "relocatable" path. Its contents will be available regardless of
whether the original path still exists or not. The contents of the path is
stored within the returned `Path` object and the folder structure is recreated
as a scratchspace if the original does not exist anymore. Calling any path
manipulation functions, such as `joinpath`, on the `Path` will return a path
to a valid folder.

`ignore` can be used to skip reading in matching files. It can be a `Regex`,
`Vector{Regex}`, or a single argument `Function` that takes the path and returns
`true` if it should be ignored, `false` otherwise.
"""
macro path(expr, ignore=:nothing)
    file = string(__source__.file)
    dir = safe_isfile(file) ? dirname(file) : pwd()
    return :($(Path)($__module__, $dir, $(esc(expr)), $(esc(ignore))))
end

struct Path <: AbstractString
    is_dir::Bool
    mod::Module
    path::String
    hash::String
    files::Dict{String,Vector{UInt8}}

    function Path(mod::Module, dir, path::AbstractString, ignore=nothing)
        path = isabspath(path) ? path : joinpath(dir, path)
        path = normpath(path)
        safe_ispath(path) || throw(ArgumentError("not a path: `$path`"))
        is_dir = isdir(path)
        dir = is_dir ? path : dirname(path)
        files = Dict{String,Vector{UInt8}}()
        ctx = SHA.SHA1_CTX()
        for (root, _, fs) in walkdir(dir), f in fs
            fullpath = joinpath(root, f)
            rel = relpath(fullpath, dir)
            if !should_ignore(rel, ignore)
                if is_dir || path == fullpath
                    include_dependency(fullpath)
                    SHA.update!(ctx, codeunits(fullpath))
                    content = read(fullpath)
                    SHA.update!(ctx, content)
                    files[rel] = content
                end
            end
        end
        return new(is_dir, mod, dir, string(Base.SHA1(SHA.digest!(ctx))), files)
    end
end

should_ignore(path::AbstractString, ::Nothing) = false
should_ignore(path::AbstractString, matcher::Function) = matcher(path)
should_ignore(path::AbstractString, re::Regex) = match(re, path) !== nothing
should_ignore(path::AbstractString, res::Vector{Regex}) = any((should_ignore(path, re) for re in res))
should_ignore(path::AbstractString, other) = error("unsupported argument type given for ignore.")

Base.show(io::IO, path::Path) = print(io, repr(getroot(path)))
Base.ncodeunits(f::Path) = ncodeunits(getpath(f))
Base.codeunit(f::Path) = codeunit(getpath(f))
Base.codeunit(f::Path, i::Integer) = codeunit(getpath(f), i)
Base.isvalid(f::Path, index::Integer) = isvalid(getpath(f), index)
Base.iterate(f::Path) = iterate(getpath(f))
Base.iterate(f::Path, state::Integer) = iterate(getpath(f), state)
Base.String(f::Path) = String(getpath(f))

function getpath(f::Path)
    if safe_ispath(f.path)
        root = getroot(f)
        # Confirm whether what's actually been returned as the root path is
        # valid, when it isn't then provide the relocated path instead.
        if safe_ispath(root)
            return root
        end
    end
    dir = Scratch.get_scratch!(f.mod, f.hash * "_" * string(hash(f.path), base=62))
    if !isempty(f.files) && !safe_ispath(joinpath(dir, first(keys(f.files))))
        cd(dir) do
            for (file, blob) in f.files
                mkpath(dirname(file))
                write(file, blob)
            end
        end
    end
    return getroot(f, dir)
end

"""
    getroot(p::Path, root = p.path)

Return the path for a `Path` object `p`, as if it were located in folder `root`.
If `p` corresponds to a folder, return `root`.
If `p` corresponds to a file `filename`, return `joinpath(root, filename)`.

!!! warning
    The result of `getroot` is not necessarily an existing path.
    This function is mostly useful to get information about `p` (e.g., its extension)
    without interacting with the filesystem.
"""
getroot(p::Path, root=p.path) = p.is_dir ? root : joinpath(root, first(keys(p.files)))

end # module
