# BoundaryIntegralEquations.jl

| **Documentation**   |  **Tests**     | **CodeCov** | **Lifecycle** |
|:--------:|:---------------:|:-------:|:-------:
|[![](https://img.shields.io/badge/docs-online-blue.svg)](https://mipals.github.io/BoundaryIntegralEquations.jl/dev/)| [![Build Status](https://github.com/mipals/BoundaryIntegralEquations.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/mipals/BoundaryIntegralEquations.jl/actions/workflows/CI.yml?query=branch%3Amain) | [![Coverage](https://codecov.io/gh/mipals/BoundaryIntegralEquations.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/mipals/BoundaryIntegralEquations.jl)| ![](https://img.shields.io/badge/Lifecycle-Unstable-yellow)| 


ViscothermalBoundaryElementMethod.jl provides the basic building blocks required for solving the Boundary Element Methods (BEM) with viscous and thermal losses in the bulk and acoustic boundary layers. Currently, it uses the collocation method for discretization of the Kirchhoff–Helmholtz integral equation found in acoustical applications

$$
c(\mathbf{y})p(\mathbf{y}) + \int_\Gamma\frac{\partial G(\mathbf{x}, \mathbf{y})}{\partial \mathbf{n} }p(\mathbf{x})\ \mathrm{d}\Gamma_\mathbf{x} = \mathrm{i}\rho ck\int_\Gamma G(\mathbf{x},\mathbf{y})v_s(\mathbf{x})\ \mathrm{d}\Gamma_\mathbf{x}.
$$

For the Fast Multipole Method this package utilizes the Julia interfaces for the Flatiron Institute Fast Multipole Libraries: [2D](https://github.com/mipals/FMM2D.jl), [3D](https://github.com/flatironinstitute/FMM3D/tree/master/julia).

For H-matrices the package utilizes the [HMatrices.jl](https://github.com/WaveProp/HMatrices.jl).

**N.B. The package is still under heavy development.**

## Installation
The package can be downloaded directly from GitHub 

```julia
using Pkg
Pkg.add(url="https://github.com/SimPr/ViscothermalBoundaryElementMethod.jl")
```

## Element types
* Continuous (Constant, Linear and Quadratic) Line elements
* Continuous (Constant, Linear and Quadratic) Triangular Elements
* Continuous (Constant, Linear and Quadratic) Quadrilateral Elements

## Mesh formats
* COMSOLs *.mphtxt* files (best for applying boundary conditions)
* .obj, .ply, .stl, .off, .2DM through [MeshIO.jl](https://github.com/JuliaIO/MeshIO.jl).

## Similar Packages
* [BEAST.jl](https://github.com/krcools/BEAST.jl): Boundary Element Analysis and Simulation Toolkit. A general toolkit, with a focus on electromagnetics. Limitations with respect to element orders and only supplies Galerkin assembly. 
* [NESSie.jl](https://github.com/tkemmer/NESSie.jl): Nonlocal Electrostatics in Structured Solvents. A specialized package written specifically for Nonlocal protein eletrostatics. 
