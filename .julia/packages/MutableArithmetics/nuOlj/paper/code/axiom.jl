function operate!!(op, args...)
    T = typeof.(args)
    if mutability(T[1], op, T...) isa IsMutable
        return operate!(op, args...)
    else
        return op(args...)
    end
end
function operate_to!!(output, op, args...)
    O = typeof(output)
    T = typeof.(args)
    if mutability(O, op, T...) isa IsMutable
        return operate_to!(output, op, args...)
    else
        return op(args...)
    end
end
