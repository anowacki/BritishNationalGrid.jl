# BritishNationalGrid

[![Build Status](https://github.com/anowacki/BritishNationalGrid.jl/workflows/CI/badge.svg)](https://github.com/anowacki/BritishNationalGrid.jl/actions)
[![Coverage Status](https://codecov.io/gh/anowacki/BritishNationalGrid.jl/graph/badge.svg?token=I27Dz8a1ER)](https://codecov.io/gh/anowacki/BritishNationalGrid.jl)

## Convert between WGS84 coordinates and British National Grid references

`BritishNationalGrid` provides the type `BNGPoint` to convert between
longitude-latitude and grid references in the [British National Grid system](https://en.wikipedia.org/wiki/Ordnance_Survey_National_Grid).
It assumes your points are geodetic longitude and latitude in decimal
degrees using the WGS84 ellipsoid.

This package is reliable to within a metre or so.  Advanced users needing
greater accuracy will probably already know how to convert between different
systems, but any additions to the package that remain easy to use will
be welcome.

## Install
```julia
julia> import Pkg

julia> Pkg.add("BritishNationalGrid")
```

## Use
Construct points in the grid using `BNGPoint`.

```julia
julia> using BritishNationalGrid

julia> p1 = BNGPoint(42513, 100231) # Full grid reference
BNGPoint{Int64}(42513, 100231)

julia> lonlat(p1) # Convert to longitude-latitude (°)
(-7.063648859478239, 50.69155306935914)

julia> p2 = BNGPoint(lon=0.32, lat=51.0) # Construct point from lon-lat
BNGPoint{Float64}(562885.4557430055, 124851.2191743746)

julia> p3 = BNGPoint(00123, 51422, "TA") # Construct from 100 km square name
BNGPoint{Int64}(500123, 451422)
```

Get a formatted grid reference:

```julia
julia> gridref(p1, 10) # 10-figure grid reference
"04251 10023"

julia> gridref(p2, 6, true) # 6-figure reference within the 100 km square TQ
"TQ 628 248"
```

Find the 100 km square in which a point lies:

```julia
julia> square(p3)
"TA"
```

### Direct (unchecked) lon-lat ↔ easting-northing conversion
The `lonlat2bng` and `bng2lonlat` functions convert directly between the
two coordinate systems, but **do not check that points are actually in the
National Grid**; therefore use these functions with caution.

```julia
julia> lon, lat = -3.183503, 55.954983;

julia> lonlat2bng(lon, lat)
(326200.06052230217, 674183.9198954724)

julia> bng2lonlat(175154, 225430)
(-5.268407238735794, 51.88199105278528)
```

## Todo
- Tie the BNGPoint type into [Geodesy.jl](https://github.com/JuliaGeo/Geodesy.jl).

## Other ways to convert to the British National Grid

- Use the Ordnance Survey's [online converter](https://www.ordnancesurvey.co.uk/gps/transformation/).  This also
  includes links to the OS's Pascal programs to do coordinate transforms.
- Use the British Geological Survey's [online converter](http://www.bgs.ac.uk/data/webservices/convertform.cfm), which also
  assumes WGS84 longitude and latitude.
