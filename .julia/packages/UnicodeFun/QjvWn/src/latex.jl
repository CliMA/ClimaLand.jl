function print_modifier(io, mod, substring)
    if mod == "^"
        to_superscript(io, substring)
    elseif mod == "_"
        to_subscript(io, substring)
    elseif mod == "bb"
        to_blackboardbold(io, substring)
    elseif mod == "bf"
        to_boldface(io, substring)
    elseif mod == "it"
        to_italic(io, substring)
    elseif mod == "cal"
        to_caligraphic(io, substring)
    elseif mod == "frak"
        to_frakture(io, substring)
    elseif mod == "mono"
        to_mono(io, substring) # leave unmodified for now
    else
        error("Modifier $mod not supported")
    end
end
"""
Base findnext doesn't handle utf8 strings correctly
"""
function utf8_findnext(A::AbstractString, v::Char, idx::Integer)
    while true
        lastidx = idx
        elem_idx = iterate(A, idx)
        elem_idx === nothing && break
        elem, idx = elem_idx
        elem == v && return lastidx
    end
    0
end

function to_latex(text)
    io = IOBuffer()
    charidx = iterate(text)
    charidx === nothing && return ""
    char, idx = charidx
    started = true
    while true
        started || (charidx = iterate(text, idx))
        started = false
        charidx === nothing && break
        char, idx = charidx
        if char in ('^', '_', '\\')
            mod = string(char)
            if mod == "\\"
                ss = SubString(text, idx, lastindex(text))
                for mod_candidate in ("bb", "bf", "it", "cal", "frak", "mono")
                    if startswith(ss, mod_candidate)
                        mod = mod_candidate
                        break
                    end
                end
                if mod == "\\" # no match was found
                    # is this a latex symbol?
                    for (k, v) in latex_symbol_map
                        if startswith(ss, k)
                            print(io, v) # replace
                            for i in 1:length(k) # move forward
                                idx = nextind(text, idx)
                            end
                            break
                        end
                    end
                    continue # ignore '\' mod
                else
                    for i in 1:length(mod) # move forward
                        idx = nextind(text, idx)
                    end
                end
            end
            char, idx = iterate(text, idx)
            if char == '{'
                i = utf8_findnext(text, '}', idx)
                if i == 0
                    error("Invalid latex. Couldn't find matching } in $(text[idx:end])")
                end
                print_modifier(io, mod, SubString(text, idx, prevind(text, i)))
                char, idx = iterate(text, i)
            else
                print_modifier(io, mod, char)
            end
        else
            print(io, char)
        end
    end
    return String(take!(io))
end
