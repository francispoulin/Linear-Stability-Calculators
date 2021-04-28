"""

Linear-Stability-Calculators
============================
Repo:      Linear-Stability-Calculators
Code:      QuasiHydrostatic/Julia/inertial_instability.jl
Model:     Quasi-Hydrostatic Model
Geometry:  Cartesian
Structure: Barotropic Bickley Jet

Nondimensional Parameters
=========================
Ro  = Rossby Number
Bu  = Burger Number
γ   = Non-traditionality parameter
Ek  = Vertical Ekman Number

Created April 22, 2021
By Francis J. Poulin

"""

abstract type AbstractGeometry end

struct Cartesian <: AbstractGeometry end
struct Spherical <: AbstractGeometry end

include("Grid.jl")                           # nothing to do here
include("Physics.jl")                        # added fₙ for non-traditional term
include("Basic_States.jl")
include("Plotting_Functions.jl")
include("Compute_Eigenvalues.jl")
