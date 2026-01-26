# PolygonOps

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://juliageometry.github.io/PolygonOps.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://juliageometry.github.io/PolygonOps.jl/dev)

The objective of this package is to provide a set of generic polygon operations. There are two assumptions: 
 - the container and points are both vector-like.
 - the first and last elements are equal. 
 
This should provide generic implementations portable and extensible for packages such as GeometryTypes, GeometryBasics, and JuliaGeo.

Please see the [docs](https://juliageometry.github.io/PolygonOps.jl/stable) for a list of exported functions.

Implements:
  - point-in-polygon
  - signed area (shoelace formula)
  - centroid

Planned:
  - simplification
  - booleans (ported from PolygonClipping.jl)
  - offset/insetting
