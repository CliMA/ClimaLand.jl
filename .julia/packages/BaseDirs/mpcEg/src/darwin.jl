function reload()
    appsupport = "~/Library/Application Support"
    # Base directories
    @setxdg DATA_HOME appsupport
    @setxdgs DATA_DIRS ["/Library/Application Support"]
    @setxdg CONFIG_HOME appsupport
    @setxdgs CONFIG_DIRS ["/Library/Application Support"]
    @setxdg STATE_HOME appsupport
    @setxdg CACHE_HOME "~/Library/Caches"
    @setxdg RUNTIME_DIR appsupport
    @setxdg BIN_HOME let
        # There's no official convention for user-specific bin files on macos,
        # so let's check the two seemingly commonly used alternatives, then fall
        # back on the system bin folder.
        usrbin = split(get(ENV, "PATH", ""), ':') ∩ [expanduser("~/.local/bin"), "/opt/local/bin"]
        if !isempty(usrbin)
            first(usrbin)
        else
            "/usr/local/bin"
        end
    end
    # User directories
    @setxdg DESKTOP_DIR "~/Desktop"
    @setxdg DOWNLOAD_DIR "~/Downloads"
    @setxdg DOCUMENTS_DIR "~/Documents"
    @setxdg MUSIC_DIR "~/Music"
    @setxdg PICTURES_DIR "~/Pictures"
    @setxdg VIDEOS_DIR "~/Videos"
    @setxdg TEMPLATES_DIR "~/Templates"
    @setxdg PUBLICSHARE_DIR "~/Public"
    # Other directories
    global FONTS_DIRS = [
        expanduser("~/Library/Fonts"),
        "/Library/Fonts",
        "/System/Library/Fonts",
        "/System/Library/Fonts/Supplemental",
        "/Network/Library/Fonts"]
    global APPLICATIONS_DIRS = ["/Applications"]
    nothing
end

"""
    simpleascii(s::String) -> String

Convert a string to a certain 'simple' ASCII form.

This is done by replacing all non-alphanumeric characters (besides `_` and `-`)
with `-`.

While this could easily be done with a regex, unlike regex
replacement this can be precompiled, reducing TTFX.
"""
Base.@assume_effects :foldable function simpleascii(s::String)
    out = Vector{UInt8}()
    sizehint!(out, ncodeunits(s))
    for b in codeunits(s)
        if UInt8('0') <= b <= UInt8('9') ||
           UInt8('A') <= b <= UInt8('Z') ||
           UInt8('a') <= b <= UInt8('z') ||
           b in (UInt8('_'), UInt8('-'))
            push!(out, b)
        else
            push!(out, UInt8('-'))
        end
    end
    String(out)
end

projectpath(p::Project, _) = projectpath(p)
projectpath(p::Project) =
    string(join(split(p.qualifier, '.') |> reverse, '.'), '.',
           simpleascii(p.org), '.',
           simpleascii(p.name), '/')
