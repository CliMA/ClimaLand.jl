module TestSyncScopes

using UnsafeAtomics: UnsafeAtomics, SyncScope
using Test

const SYNCSCOPES = [:none, :singlethread]

function check_syncscope(name::Symbol)
    sync = getfield(UnsafeAtomics, name)
    @test sync isa SyncScope
    @test string(sync) == sprint(print, sync) == string(sync)
    @test repr(sync) == "UnsafeAtomics.$name"
end

function test_syncscope()
    @testset for name in SYNCSCOPES
        check_syncscope(name)
    end
end

end  # module