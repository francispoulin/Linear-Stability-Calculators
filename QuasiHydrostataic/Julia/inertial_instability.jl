#!/usr/bin/env julia --startup-file=no

"""
Linear-Stability-Calculators: Quasi-Hydrostatic Model
=====================================================

Linear Stability Calculation of an inertially unstable jet using spectral methods

The basic state is that of a barotropic Bickley jet (to be generalized).
The perturbation is assumed to be periodic in both the zonal (x) and vertical (z) directions.

Calculations are done in serial.  See inertial_instability_mpi.py for parallel execution.
"""

using Revise

abstract type AbstractMesh end

struct Cheb <: AbstractMesh end
struct FiniteDifferenceOrder2  <: AbstractMesh end

include("src/parameters.jl")
include("src/Compute_Eigenvalues.jl")
include("src/Plotting_Functions.jl")

files   = Files()
grid    = Grid(Ly = 1000e3, Lz = 3e3, Ny = 10, θ₀ = π/32, method=Cheb())
physics = Physics(N = 1e-2, ν = 0.26, θ₀ = grid.θ₀, NT = 1)
jet     = Jet(grid, physics)

Output_Parameters(grid, physics, jet, files.json)

### Number of eigenvalues to store and wavenumbers
Neigs  = 1                            
dk, Nk = 1e-6, 41
dm, Nm = 1e-4, 150

ks = collect(0 :dk:(Nk-1)*dk)
ms = collect(dm:dm:Nm*dm)

ωs    = zeros(Complex, (Nk, Nm, Neigs))
modes = zeros(Complex, (Nk, Nm, 3*grid.Ny+1, Neigs))

for (ik, k) in enumerate(ks)
    print("ik = ", ik, "\n")
    for (im, m) in enumerate(ms)
        A = build_matrix(k, m, grid, physics, jet)

        ωs[ik, im, :], modes[ik, im, :, :] = compute_spectrum(A, Neigs)
    end
end

save_spectrum(ωs, modes, ks, ms, Neigs, grid.y, grid.Ny, files.nc)

plot_growth_slice(files.nc, files.json, files.plotslicem)
plot_growth(      files.nc, files.json, files.plotgrowth)

plot_modes_1D(files.nc, files.json, files.plotmodes1D, Neigs)
plot_modes_2D(files.nc, files.json, files.plotmodes2D, Neigs)
