# ChebPade.jl

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://mjp98.github.io/ChebPade.jl/dev)
[![Build Status](https://github.com/mjp98/ChebPade.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/mjp98/ChebPade.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/mjp98/ChebPade.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/mjp98/ChebPade.jl)

A Julia package to compute the Chebyshev-Pade approximation from a Chebyshev polynomial expansion, based on the methods implemented in [Chebfun](https://github.com/chebfun/chebfun) [Chebfun license](https://github.com/chebfun/chebfun/blob/master/LICENSE.txt).

## Installation

This package is unregistered, so may be added by 

```julia
(v1.7) pkg> add https://github.com/mjp98/ChebPade.jl
```

This package supports Julia v1.7 and later.

## Usage

```julia
julia> using ChebPade, ApproxFun
```
### Example

Construct a Chebyshev-Pade approximant of the exponential function of type (2,2), returning the numerator and denominator of type ApproxFun.Fun

```julia
julia> f = Fun(exp,-1..1)
Fun(Chebyshev(-1..1),[1.2660658777520084, 1.1303182079849703, 0.27149533953407656, 0.04433684984866379, 0.0054742404420936785, 0.0005429263119139232, 4.497732295427654e-5, 3.19843646253781e-6, 1.992124804817033e-7, 1.1036771869970875e-8, 5.505896578301994e-10, 2.4979607981699334e-11, 1.0391104209722668e-12, 3.993680386393805e-14])

julia> p,q = chebpade(f,2,2)
(Fun(Chebyshev(-1..1),[1.0020016925359065, 0.48492995699568886, 0.04024039039994624]), Fun(Chebyshev(-1..1),[1.0, -0.47643278089764673, 0.038277919301135936]))

julia> p,q = chebpade(f,2,2)
(Fun(Chebyshev(-1..1),[1.0020016925359065, 0.48492995699568886, 0.04024039039994624]), Fun(Chebyshev(-1..1),[1.0, -0.47643278089764673, 0.038277919301135936]))

julia> using LinearAlgebra; p./q - f |> norm
9.313144965137737e-5
```

## To do:

 - improve matrix construction to reduce allocations
 - optimise matrix inversion, e.g. inplace
