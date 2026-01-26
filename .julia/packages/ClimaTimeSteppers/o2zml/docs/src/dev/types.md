# Types

In this section, we print out the type hierarchies of some classes via code snippets.

## Algorithms

```@example
import AbstractTrees as AT
import InteractiveUtils as IU
import ClimaTimeSteppers as CTS
AT.children(x::Type) = IU.subtypes(x)
AT.print_tree(CTS.DistributedODEAlgorithm)
```

## Tableaus

```@example
import AbstractTrees as AT
import InteractiveUtils as IU
import ClimaTimeSteppers as CTS
AT.children(x::Type) = IU.subtypes(x)
AT.print_tree(CTS.AbstractAlgorithmName)
```
