function emit_method_check(mod::Module, name::Symbol)
    msg = "The object $name is not owned by $mod, you cannot declare it as public."
    return quote
        if $mod !== $Base.parentmodule($name)
            $Base.error($msg)
        end
    end
end
