macro declare_file_product(product_name)
    get_path_name = Symbol(string("get_", product_name, "_path"))
    path_name = Symbol(string(product_name, "_path"))
    return esc(quote
        # These will be filled in by init_file_product().
        @static if $global_typeassert_available
            $(product_name)::String = ""
            $(path_name)::Union{Nothing,String} = $(emit_preference_path_load(string(product_name, "_path")))
        else
            $(product_name) = ""
            $(path_name) = $(emit_preference_path_load(string(product_name, "_path")))
        end
        function $(get_path_name)()
            return $(path_name)::String
        end
    end)
end

macro init_file_product(product_name, product_path)
    path_name = Symbol(string(product_name, "_path"))
    return esc(quote
        global $(path_name)
        # FileProducts are very simple, and we maintain the `_path` suffix version for consistency
        if $(path_name) === nothing
            $(path_name) = joinpath(artifact_dir, $(product_path))
        end
        global $(product_name) = $(path_name)
    end)
end
