```@meta
CurrentModule = ADTypes
CollapsedDocStrings = true
```

# ADTypes.jl

Documentation for [ADTypes.jl](https://github.com/SciML/ADTypes.jl/).

```@docs
ADTypes
AbstractADType
```

## Dense AD

### Forward mode

Algorithmic differentiation:

```@docs
AutoForwardDiff
AutoPolyesterForwardDiff
AutoMooncakeForward
```

Finite differences:

```@docs
AutoFiniteDiff
AutoFiniteDifferences
```

Taylor mode:

```@docs
AutoGTPSA
AutoTaylorDiff
```

### Reverse mode

```@docs
AutoMooncake
AutoReverseDiff
AutoTracker
AutoZygote
```

### Forward or reverse mode

```@docs
AutoEnzyme
AutoChainRules
AutoDiffractor
```

### Symbolic mode

```@docs
AutoFastDifferentiation
AutoSymbolics
```

## Sparse AD

```@docs
AutoSparse
ADTypes.dense_ad
```

### Sparsity detector

```@docs
ADTypes.sparsity_detector
ADTypes.AbstractSparsityDetector
ADTypes.jacobian_sparsity
ADTypes.hessian_sparsity
ADTypes.NoSparsityDetector
ADTypes.KnownJacobianSparsityDetector
ADTypes.KnownHessianSparsityDetector
```

### Coloring algorithm

```@docs
ADTypes.coloring_algorithm
ADTypes.AbstractColoringAlgorithm
ADTypes.column_coloring
ADTypes.row_coloring
ADTypes.symmetric_coloring
ADTypes.NoColoringAlgorithm
```

## No automatic differentiation

```@docs
NoAutoDiff
NoAutoDiffSelectedError
```

## Modes

```@docs
ADTypes.mode
ADTypes.AbstractMode
ADTypes.ForwardMode
ADTypes.ForwardOrReverseMode
ADTypes.ReverseMode
ADTypes.SymbolicMode
```

## Miscellaneous

```@docs
ADTypes.Auto
```

## Deprecated

```@docs
AutoTapir
```
