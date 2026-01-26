import CommonDataModel as CDM

struct DummyEmptyDataset <: CDM.AbstractDataset
end


dd = DummyEmptyDataset();
@test length(CDM.dimnames(dd)) == 0
@test_throws Exception CDM.dim(dd,"does_not_exist")
@test_throws Exception CDM.defDim(dd,"does_not_exist",1)

@test length(CDM.attribnames(dd)) == 0
@test_throws Exception CDM.attrib(dd,"does_not_exist")
@test_throws Exception CDM.defAttrib(dd,"does_not_exist",1)

@test length(CDM.varnames(dd)) == 0
@test_throws Exception CDM.variable(dd,"does_not_exist")
@test_throws Exception CDM.defVar(dd,"does_not_exist",Int32,())

@test CDM.path(dd) == ""

@test length(CDM.groupnames(dd)) == 0
# not available in julia 1.6
#@test_throws "no group" CDM.group(dd,"does_not_exist")
#@test_throws "unimplemented" CDM.defGroup(dd,"does_not_exist")
@test_throws Exception CDM.group(dd,"does_not_exist")
@test_throws Exception CDM.defGroup(dd,"does_not_exist")
