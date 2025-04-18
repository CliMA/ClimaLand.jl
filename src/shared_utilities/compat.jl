import ClimaCore

@static if pkgversion(ClimaCore) < v"0.14.30"
    function compat_set_mask!(space)
        return nothing
    end
else
    function compat_set_mask!(space)
        ClimaCore.Spaces.set_mask!(space, ClimaCore.Fields.ones(space))
    end
end

@static if pkgversion(ClimaCore) < v"0.14.30"
    function compat_add_mask()
        return (;)
    end
else
    function compat_add_mask()
        (; enable_mask = true)
    end
end
