abstract type AbstractGeometry end

struct Cartesian <: AbstractGeometry end
struct Spherical <: AbstractGeometry end

include("cheb.jl")
include("parameters.jl")

#=
"""
Cartesian version
"""
function test_dispatch(geometry::Cartesian)
    print("Cartesian")
end

"""
Spherical version
"""
function test_dispatch(geometry::Spherical)
    print("Spherical")
end
=#