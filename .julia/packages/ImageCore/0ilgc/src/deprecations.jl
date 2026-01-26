const warned_once = Ref(false)
function forced_depwarn(msg, sym)
    opts = Base.JLOptions()
    if !warned_once[] && !(opts.depwarn == 1)
        @warn msg
        @info """It is recommended that you fix this now to avoid breakage when a new version is released and this warning is removed.
                    Tip: to see all deprecation warnings together with code locations, launch Julia with `--depwarn=yes` and rerun your code."""
        warned_once[] = true
    else
        Base.depwarn(msg, sym)
    end
    return nothing
end

# a perhaps "permanent" deprecation
Base.@deprecate_binding permuteddimsview PermutedDimsArray
