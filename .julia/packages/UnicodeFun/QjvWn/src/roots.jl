function to_root(io::IO, root, val)
    if root == 1
        write(io, string(val))
        return
    end
    if root == 3
        print(io, '∛')
    elseif root == 4
        print(io, '∜')
    else
        root != 2 && to_superscript(io, root)
        print(io, '√')
    end
    to_overline(io, string(val))
end
to_root(io::IO, val) = to_root(io, 2, val)

function to_root(root, val)
    sprint() do io
        to_root(io, root, val)
    end
end
to_root(val) = to_root(2, val)

