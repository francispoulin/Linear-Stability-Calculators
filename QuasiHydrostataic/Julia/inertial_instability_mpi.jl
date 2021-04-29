#!/usr/bin/env julia --startup-file=no

"""
Linear-Stability-Calculators: Quasi-Hydrostatic Model
=====================================================

Linear Stability Calculation of an inertially unstable jet using spectral methods

The basic state is that of a barotropic Bickley jet (to be generalized).
The perturbation is assumed to be periodic in both the zonal (x) and vertical (z) directions.

Calculations are done in parallel.  See inertial_instability.py for serial execution.
"""

using MPI
using Printf

MPI.Init()

comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)
comm_size = MPI.Comm_size(comm)

abstract type AbstractMesh end

struct Cheb <: AbstractMesh end
struct FiniteDifferenceOrder2  <: AbstractMesh end

include("src/parameters.jl")
include("src/mpi_stuff.jl")
include("src/Compute_Eigenvalues.jl")
#include("src/Plotting_Functions.jl")

files   = Files()
grid    = Grid(Ly = 1000e3, Lz = 3e3, Ny = 10, θ₀ = π/32, method=Cheb())
physics = Physics(N = 1e-2, ν = 0.26, θ₀ = grid.θ₀, NT = 1)
jet     = Jet(grid, physics)

rank == 0 ? Output_Parameters(grid, physics, jet, files.json) : nothing

### Number of eigenvalues to store and wavenumbers
Neigs  = 10                            
dk, Nk = 1e-6, 11 #41
dm, Nm = 1e-4, 20 #150

rank == 0 ? ks = Float64[(ik-1)*dk for ik in 1:grid.Ny+1] : k = nothing
#rank == 0 ? @printf("original ks = %s \n", ks) : nothing

ms = collect(dm:dm:Nm*dm)

Nk_local, iks_ends = array_split(Nk, comm_size, comm)

ks_local = zeros(Float64, Nk_local[rank+1])
MPI.Scatterv!(rank == 0 ? VBuffer(ks, Nk_local, iks_ends) : nothing, ks_local, 0, comm)
#@printf("rank = %s  k_local = %s \n", rank, ks_local)

ωs_local    = zeros(Complex{Float64}, (Nk_local[rank+1], Nm, Neigs))
modes_local = zeros(Complex{Float64}, (Nk_local[rank+1], Nm, 3*grid.Ny+1, Neigs))

for (ik, k) in enumerate(ks_local)
    print("ik = ", ik, "\n")
    for (im, m) in enumerate(ms)
        A = build_matrix(k, m, grid, physics, jet)

        ωs_local[ik, im, :], modes_local[ik, im, :, :] = compute_spectrum(A, Neigs)
    end
end

### Gatter array
rank == 0 ? ωs = zeros(Complex{Float64}, (Nk, Nm, Neigs)) : nothing

rank == 0 ? @printf("rank = %s and size ωs =%s\n", rank, size(ωs)) : nothing
@printf("rank = %s, size ωs_local = %s\n", rank, size(ωs_local))

MPI.Gatherv!(ωs_local, rank == 0 ? VBuffer(ωs, Nk_local, iks_ends) : nothing, 0, comm)
#rank == 0 ? @printf("gathered k = %s\n", ωs) : nothing

#=
-> must gather arrays then save as normal and plot as normal 
save_spectrum(ωs, modes, ks, ms, Neigs, grid.y, grid.Ny, files.nc)

plot_growth_slice(files.nc, files.json, files.plotslicem)
plot_growth(      files.nc, files.json, files.plotgrowth)

plot_modes_1D(files.nc, files.json, files.plotmodes1D, Neigs)
plot_modes_2D(files.nc, files.json, files.plotmodes2D, Neigs)
=#
