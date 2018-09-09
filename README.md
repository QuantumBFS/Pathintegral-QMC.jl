# Pathintegral-QMC.jl
Path integral Monte Carlo for solving Transverse Ising Model et. al.

The initial commit is a transription of Stephan Humeniuk's excellent Fortran code. It applies SSE on transverse field Ising model (TFI) to calculate ground state energy and magnetic properties.

We start from simplest path integral Monte Carlo, but the framework is designed for extensibility.

## How to use
(under development, think before using it)
```
julia main.jl
```
Julia version should be 1.0.

## Authors
* Stephan Humeniuk
* JinGuo Liu
