function emit_stub(mod::Module)
    if !isdefined(mod, JIEKO_STUB)
        return :(const $JIEKO_STUB = $JiekoStub())
    end
    return nothing
end
