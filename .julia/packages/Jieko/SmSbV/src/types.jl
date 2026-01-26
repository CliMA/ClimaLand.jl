abstract type Captured end

struct CapturedStruct <: Captured
    mod::Module
    name::Symbol
    doc::String
end

struct CapturedFunction <: Captured
    mod::Module
    name::Symbol
    sig::Type
    rettype::Type
    doc::String
end

struct CapturedMacro <: Captured
    mod::Module
    name::Symbol
    doc::String
end

struct CapturedConst <: Captured
    mod::Module
    name::Symbol
    doc::String
end

struct JiekoStub
    consts::Dict{Symbol,CapturedConst}
    macros::Dict{Symbol,CapturedMacro}
    structs::Dict{Symbol,CapturedStruct}
    interface::Dict{Type,CapturedFunction}
end

function JiekoStub()
    return JiekoStub(
        Dict{Symbol,CapturedConst}(),
        Dict{Symbol,CapturedMacro}(),
        Dict{Symbol,CapturedStruct}(),
        Dict{Type,CapturedFunction}()
    )
end

function allcaptured(stub::JiekoStub)
    return Iterators.flatten([
        values(stub.consts),
        values(stub.macros),
        values(stub.structs),
        values(stub.interface)
    ])
end

function Base.setindex!(stub::JiekoStub, value::CapturedConst, key::Symbol)
    stub.consts[key] = value
    return stub
end

function Base.setindex!(stub::JiekoStub, value::CapturedMacro, key::Symbol)
    stub.macros[key] = value
    return stub
end

function Base.setindex!(stub::JiekoStub, value::CapturedStruct, key::Symbol)
    stub.structs[key] = value
    return stub
end

function Base.setindex!(stub::JiekoStub, value::CapturedFunction, key::Type)
    stub.interface[key] = value
    return stub
end
