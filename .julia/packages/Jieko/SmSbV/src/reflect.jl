"""
$DEF

Return the interface stub storage for a module.
"""
@pub function stub(mod::Module)::JiekoStub
    if isdefined(mod, JIEKO_STUB)
        return getfield(mod, JIEKO_STUB)
    else
        return JiekoStub()
    end
end
