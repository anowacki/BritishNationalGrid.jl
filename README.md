# BritishNationalGrid

[![Build Status](https://img.shields.io/travis/anowacki/BritishNationalGrid.jl.svg?style=flat-square&label=linux)](https://travis-ci.org/anowacki/BritishNationalGrid.jl)
[![Coverage Status](https://coveralls.io/repos/github/anowacki/BritishNationalGrid.jl/badge.svg?branch=master)](https://coveralls.io/github/anowacki/BritishNationalGrid.jl?branch=master)

## Convert between WGS84 coordinates and British National Grid references

`BritishNationalGrid` provides the type `BNGPoint` to convert between
longitude-latitude and grid references.

## Install
```julia
julia> Pkg.clone("https://github.com/anowacki/BritishNationalGrid.jl")
```

**NB**: Automatic installation of the Proj4 library [doesn't seem to be working
on Windows yet](https://github.com/JuliaGeo/Proj4.jl/issues/8)

## Use
Construct points in the grid using `BNGPoint`.

```julia
julia> using BritishNationalGrid

julia> p1 = BNGPoint(42513, 100231) # Full grid reference
BritishNationalGrid.BNGPoint{Int64}(42513, 100231)

julia> lonlat(p1) # Convert to longitude-latitude (°)
(-7.063648859478239, 50.691553069358555)

julia> p2 = BNGPoint(lon=0.32, lat=51.0) # Construct point from lon-lat
BritishNationalGrid.BNGPoint{Float64}(562885.4557430055, 124851.2191743746)

julia> p3 = BNGPoint(00123, 51422, "TA") # Construct from 100 km square name
BritishNationalGrid.BNGPoint{Int64}(500123, 451422)
```

Get a formatted grid reference:

```julia
julia> gridref(p1, 10) # 10-figure grid reference
"04251 10023"

julia> gridref(p2, 6, true) # 8-figure reference within the 100 km square
"TQ 628 248"
```

Find the 100 km square in which a point lies:

```julia
julia> square(p3)
"TA"
```

## Todo
- Tie the BNGPoint type into [Geodesy.jl](https://github.com/JuliaGeo/Geodesy.jl).
