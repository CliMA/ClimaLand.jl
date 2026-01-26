
## Data types

In order to implement a new dataset based `CommonDataModel.jl`
one has to create two types derived from:

1. [`AbstractVariable`](@ref CommonDataModel.AbstractVariable): a variable with named dimension and metadata
2. [`AbstractDataset`](@ref CommonDataModel.AbstractDataset): a collection of variable with named dimension, metadata and sub-groups. The sub-groups are also `AbstractDataset`.


`CommonDataModel.jl` also provides a type `CFVariable` which wraps a type derived from `AbstractVariable` and applies the scaling described in
[`cfvariable`](@ref CommonDataModel.cfvariable).

Overview of methods:

|            | get names                                     | get values                              | set value                        | property |
|------------|-----------------------------------------------|-----------------------------------------|-------------------------------------------|--------|
| Dimensions | [`dimnames`](@ref CommonDataModel.dimnames)       | [`dim`](@ref CommonDataModel.dim)           | [`defDim`](@ref CommonDataModel.defDim)       | `dim`    |
| Attributes | [`attribnames`](@ref CommonDataModel.attribnames) | [`attrib`](@ref CommonDataModel.attrib)     | [`defAttrib`](@ref CommonDataModel.defAttrib) | `attrib` |
| Variables  | [`varnames`](@ref CommonDataModel.varnames)    | [`variable`](@ref CommonDataModel.variable) | [`defVar`](@ref CommonDataModel.defVar)       | -      |
| Groups     | [`groupnames`](@ref CommonDataModel.groupnames)   | [`group`](@ref CommonDataModel.group)       | [`defGroup`](@ref CommonDataModel.defGroup)   | `group`  |

For read-only datasets, the methods in "set value" column are not to be implemented.
Attributes can also be delete with the [`delAttrib`](@ref CommonDataModel.delAttrib) functions.

Every struct deriving from `AbstractDataset` have automaticaly the special properties `dim`, `attrib` and `group` which act like dictionaries (unless a field with this name already exists).
For `attrib`, calls to `keys`, `getindex` and `setindex!`, `delete!` are dispated to `attribnames`, `attrib`,`defAttrib`, and `delAttrib` respectively (and likewise for other properties). For example:

``` julia
using NCDatasets
ds = NCDataset("file.nc")
# setindex!(ds.attrib,...) here automatically calls defAttrib(ds,...)
ds.attrib["title"] = "my amazing results";
```
Variables can be accessed by directly indexing the `AbstractDataset`.

Every struct deriving from `AbstractVariable` has the properties `dim`, and `attrib`.

Current functionalities of CommonDataModel include:
* virtually concatenating files along a given dimension
* create a virtual subset (([`view`](@ref Base.view))) by indices or by values of coordinate variables ([`select`](@ref CommonDataModel.select), [`@select`](@ref CommonDataModel.@select))
* group, map and reduce a variable ([`groupby`](@ref CommonDataModel.groupby), [`@groupby`](@ref CommonDataModel.@groupby), [`rolling`](@ref CommonDataModel.rolling))



## API

```@autodocs
Modules = [CommonDataModel, CommonDataModel.CatArrays]
```
