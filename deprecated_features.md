# Deprecated Features List

This file lists deprecated features and the recommended alternative practice

- ## `root_distribution`

    The `root_distribution` in `PlantHydraulicsParameters` is replaced by `rooting_depth`.
    If using a `root_distribution` function of the form P(x) = 1/d e^(z/d), then `rooting_depth`
    is equivalent to d.
