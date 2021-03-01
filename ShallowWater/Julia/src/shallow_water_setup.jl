abstract type AbstractGeometry end

struct Cartesian <: AbstractGeometry end
struct Spherical <: AbstractGeometry end

include("Grid.jl")
include("Physics.jl")
include("Basic_States.jl")
include("Plotting_Functions.jl")
include("Compute_Eigenvalues.jl")
