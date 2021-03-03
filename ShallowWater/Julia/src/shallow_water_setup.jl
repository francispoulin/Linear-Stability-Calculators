"""

Linear-Stability-Calculators
============================
Repo:      Linear-Stability-Calculators
Code:      ShallowWater/Julia/linear_stability_shallow_water.jl
Model:     Rotating Shallow Water
Geometry:  Cartesian/Spherical
Structure: Bickley Jet

Nondimensional Parameters
=========================
Fr = Froude Number
Ro = Rossby Number

Created Feb 26, 2021
By Francis J. Poulin

"""

abstract type AbstractGeometry end

struct Cartesian <: AbstractGeometry end
struct Spherical <: AbstractGeometry end

include("Grid.jl")
include("Physics.jl")
include("Basic_States.jl")
include("Plotting_Functions.jl")
include("Compute_Eigenvalues.jl")
