"""
### BritishNationalGrid

Convert between longitude and latitude and coordinates of the
British National Grid.

The `BNGPoint` type and constructor is used to create points with easting
and northings  These are given relative to the grid's false origin southwest
of the Isles of Scilly.

Exported methods:

- `BNGPoint`: Construct a new point on the grid
- `gridref`: Return a string with an n-figure grid reference
- `lonlat`: Convert a grid point to WGS84 longitude and latitude
- `square`: Determine which National Grid 100 km square a point is in

See the documentation for each method to learn more.
"""
module BritishNationalGrid

using Artifacts: @artifact_str
using Proj: Proj, Transformation
using Printf: Format, format

export
    BNGPoint,
    gridref,
    lonlat,
    square

const BNG_STRING = "EPSG:27700"
const WGS84_STRING = "EPSG:4326"

const FROM_WGS84 = Ref{Transformation}()
const FROM_BNG = Ref{Transformation}()

function __init__()
    # Allow Proj to find the grid needed to convert from OSGB to ETRS
    # by adding the artifact directory to its search path
    Proj.proj_context_set_search_paths(2, [Proj.PROJ_DATA[], artifact"ostn15_grid"])
    FROM_WGS84[] = Transformation(WGS84_STRING, BNG_STRING; always_xy=true)
    FROM_BNG[] = Transformation(BNG_STRING, WGS84_STRING; always_xy=true)
end

"""
#### Type

    BNGPoint{T<:Real}

A struct holding the easting and northing of a point within the British
National Grid.

#### Constructors

    BNGPoint(e, n)

Provide eastings `e` and northings `n` in m from the grid origin.

    BNGPoint(e, n, square)

Provide a point within a 100 km square and its name as a two-character string,
creating a point with a full reference

```jldoctest
julia> BNGPoint(101, 12345, "OV")
BritishNationalGrid.BNGPoint{Int64}(500101, 512345)

```

    BNGPoint(; lon=0, lat=0)

Convert a WGS84 longitude `lon` and latitude `lat` in degrees, into a grid point.

```jldoctest
julia> BNGPoint(lon=-1.54, lat=55.5)
BNGPoint{Float64}(429158.6966862687, 623009.0927592682)

```

#### Acessing fields
A `BNGPoint` contains the following fields:
- `e`: Easting in m
- `n`: Northing in m

These fields are part of the API, and therefore it is safe and expected
for users to extract the eastings and northings from a `BNGPoint` like so:

```jldoctest
julia> p = BNGPoint(lon=-1.54, lat=55.5)
BNGPoint{Float64}(429158.6966862687, 623009.0927592682)

julia> p.e, p.n
(429158.6966862687, 623009.0927592682)
"""
struct BNGPoint{T<:Real}
    e::T
    n::T
    function BNGPoint{T}(e, n) where {T<:Real}
        in_grid(e, n) ||
            throw(ArgumentError("Easting and northing are not within the grid. "*
                "(Are: $e, $n)"))
        new(e, n)
    end
end

BNGPoint(e::T1, n::T2) where {T1<:Real,T2<:Real} = BNGPoint{promote_type(T1, T2)}(e, n)

function BNGPoint(e::T1, n::T2, sq::String) where {T1, T2}
    0 <= e < 100_000. && 0. <= n < 100_000. ||
        throw(ArgumentError("Easting and/or northing are not within a 100 km " *
            "square.  (Are: $e, $n.)"))
    sq in SQUARE_NAMES || throw(ArgumentError("'$sq' is not a valid square name"))
    iN, iE = -1, -1
    for i in eachindex(SQUARE_NAMES)
        if SQUARE_NAMES[i] == sq
            iN, iE = Tuple(CartesianIndices(SQUARE_NAMES)[i])
            break
        end
    end
    @assert (iN, iE) != (-1, -1)
    e += (iE-1)*100_000.0
    n += (iN-1)*100_000.0
    BNGPoint{promote_type(T1, T2)}(e, n)
end
function BNGPoint(v::Union{AbstractVector{T},Tuple{T,T}}) where {T<:Real}
    length(v) == 2 || throw(ArgumentError("Length-2 vector or tuple `v` required"))
    BNGPoint(v[1], v[2])
end
BNGPoint(; lon=0.0, lat=0.0) = BNGPoint.(lonlat2bng.(lon, lat))

"""
    gridref(p::BNGPoint, n; square=false, sep=' ')

Return a string giving an `n`-figure grid reference.  By default, a full reference
is given.  If `square` is `true`, then supply the 100 km square name first, then
the reference within that square.  The square, eastings and northings are
separated by `sep`.

```jldoctest
julia> gridref(BNGPoint(429157, 623009), 8, square=true, sep="_")
"NU_2915_2300"

```
"""
function gridref(p::BNGPoint, n::Integer=8; square::Bool=false, sep=' ')
    2 <= n <= (square ? 10 : 12) ||
        throw(ArgumentError("Grid references must be given to between 2 " *
            "and 12 digits without the name of the 100 km square, or 2 to 10 " *
            "with (asked for $n)"))
    n%2 == 0 || throw(ArgumentError("Number of figures must be even"))
    n = n÷2
    east, north = p.e, p.n
    if square
        east %= 100_000
        north %= 100_000
    end
    divisor = 10.0^(max(6 - (square ? n+1 : n), 0))
    east = floor(Int, east/divisor)
    north = floor(Int, north/divisor)
    fmt = Format("%0$(n)d")
    io = IOBuffer()
    if square
        write(io, BritishNationalGrid.square(p))
        if !isempty(sep)
            write(io, sep)
        end
    end
    format(io, fmt, east)
    if !isempty(sep)
        write(io, sep)
    end
    format(io, fmt, north)
    String(take!(io))
end

"""
    lonlat(p::BNGPoint) -> lon, lat

Return the WGS84 longitude `lon` and latitude `lat` in decimal degrees for the point `p`.
"""
lonlat(p::BNGPoint) = bng2lonlat(p.e, p.n)

"""
    square(p::BNGPoint) -> XX::String

Return a two-character string `XX` containing the name of the 100 km-by 100 km
square in which is located the point `p`.

```jldoctest
julia> using BritishNationalGrid

julia> BritishNationalGrid.square(BNGPoint(200_000, 1_000_000))
"HX"

```
"""
square(p::BNGPoint) = _square(p.e, p.n)


# Internal routines
"""
    lonlat2bng(lon, lat) -> easting, northing
    lonlat2bng(lon::AbstractArray, lat::AbstractArray) -> A::Array

Transform from longitude and latitude in WGS84 into BNG easting and northing (m).
The first form does so for scalars and returns a tuple; the second form does
so for length-n arrays and returns a n-by-2 array where the first column is the
easting, and the second is the northing.
"""
function lonlat2bng(lon::Real, lat::Real)
    e, n = FROM_WGS84[](lon, lat)
    e, n
end

function bng2lonlat(e::T1, n::T2) where {T1<:Real, T2<:Real}
    lon, lat = FROM_BNG[](e, n)
    lon, lat
end

"""
    in_grid(e, n) -> ::Bool

Return `true` if the easting `e` and northing `n` (in m) are within the
British National Grid, and `false` otherwise.
"""
in_grid(e, n) = 0. <= e < 700_000. && 0. <= n < 1_300_000.

function _square(e, n)
    in_grid(e, n) || throw(ArgumentError("Point is outside the grid"))
    iE = floor(Int, e/100_000.) + 1
    iN = floor(Int, n/100_000.) + 1
    SQUARE_NAMES[iN,iE]
end

"""
    square_names() -> names::Array{String,2}

Build the two-letter codes of each 100 km-b-100 km square of the grid.
Access the `names` by:

```julia
julia> using BritishNationalGrid

julia> easting, northing = 200_000, 1_000_000
(200000, 1000000)

julia> squares = BritishNationalGrid.square_names()
13×7 Array{String,2}:
 "SV"  "SW"  "SX"  "SY"  "SZ"  "TV"  "TW"
 "SQ"  "SR"  "SS"  "ST"  "SU"  "TQ"  "TR"
 "SL"  "SM"  "SN"  "SO"  "SP"  "TL"  "TM"
 "SF"  "SG"  "SH"  "SJ"  "SK"  "TF"  "TG"
 "SA"  "SB"  "SC"  "SD"  "SE"  "TA"  "TB"
 "NV"  "NW"  "NX"  "NY"  "NZ"  "OV"  "OW"
 "NQ"  "NR"  "NS"  "NT"  "NU"  "OQ"  "OR"
 "NL"  "NM"  "NN"  "NO"  "NP"  "OL"  "OM"
 "NF"  "NG"  "NH"  "NJ"  "NK"  "OF"  "OG"
 "NA"  "NB"  "NC"  "ND"  "NE"  "OA"  "OB"
 "HV"  "HW"  "HX"  "HY"  "HZ"  "JV"  "JW"
 "HQ"  "HR"  "HS"  "HT"  "HU"  "JQ"  "JR"
 "HL"  "HM"  "HN"  "HO"  "HP"  "JL"  "JM"

julia> squares[floor(Int, northing/100_000)+1, floor(Int, easting/100_000)+1]
"HX"

```
"""
function square_names()
    names = Array{String}(undef, 13, 7)
    letter2 = ["A" "B" "C" "D" "E"
               "F" "G" "H" "J" "K"
               "L" "M" "N" "O" "P"
               "Q" "R" "S" "T" "U"
               "V" "W" "X" "Y" "Z"]
    names[1:3,1:5] .= "H".*letter2[3:5,:]
    names[4:8,1:5] .= "N".*letter2
    names[9:13,1:5] .= "S".*letter2
    names[1:3,6:7] .= "J".*letter2[3:5,1:2]
    names[4:8,6:7] .= "O".*letter2[:,1:2]
    names[9:13,6:7] .= "T".*letter2[:,1:2]
    names[end:-1:1,:]
end
const SQUARE_NAMES = square_names()


end # module
