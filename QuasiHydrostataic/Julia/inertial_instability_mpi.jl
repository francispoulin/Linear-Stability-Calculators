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
include("src/Plotting_Functions.jl")

files   = Files()
grid    = Grid(Ly = 1000e3, Lz = 3e3, Ny = 20, θ₀ = π/32, method=Cheb())
physics = Physics(N = 1e-2, ν = 0.26, θ₀ = grid.θ₀, NT = 1)
jet     = Jet(grid, physics)

rank == 0 ? Output_Parameters(grid, physics, jet, files.json) : nothing

Neigs  = 1                            
dk, Nk = 1e-6, 2 #11 #41
dm, Nm = 1e-4, 4 #150

rank == 0 ? ks = Float64[(ik-1)*dk for ik in 1:Nk] : k = nothing

ms = collect(dm:dm:Nm*dm)

Nk_local, iks_ends = array_split(Nk, comm_size, comm)

ks_local = zeros(Float64, Nk_local[rank+1])
MPI.Scatterv!(rank == 0 ? VBuffer(ks, Nk_local, iks_ends) : nothing, ks_local, 0, comm)

ωs_local    = zeros(Complex{Float64}, (Nk_local[rank+1], Nm, Neigs))
modes_local = zeros(Complex{Float64}, (Nk_local[rank+1], Nm, 3*grid.Ny+1, Neigs))

for (ik, k) in enumerate(ks_local)
    #if rank == 0
    #    @printf("ik = %s \n", ik)
    #end
    for (im, m) in enumerate(ms)
        A = build_matrix(k, m, grid, physics, jet)

        ωs_local[ik, im, :], modes_local[ik, im, :, :] = compute_spectrum(A, Neigs)

        #print("k = ", ks[ik], " m = ", ms[im], " growth = ", imag(ωs_local[ik, im, 1]), "\n")
    end
end

for (ik, k) in enumerate(ks_local)
    @printf("r=%s,  k=%s, ωs=%s \n", rank, round(k, digits=6), round.(ωs_local[ik, :, 1]/physics.fz, digits=6))
end
print("\n")

if rank == 0
    ωs    = zeros(Complex{Float64}, (Nk, Nm, Neigs)) 
    modes = zeros(Complex{Float64}, (Nk, Nm, 3*grid.Ny+1, Neigs))
end

counts = Nk_local * Nm * Neigs
displs = cumsum(append!([0], counts))[1:comm_size]
MPI.Bcast!(counts, 0, comm)
MPI.Bcast!(displs, 0, comm)
MPI.Gatherv!(ωs_local, rank == 0 ? VBuffer(ωs, counts, displs) : nothing, 0, comm)

if rank == 0
    for (ik, k) in enumerate(ks)
        @printf("r=%s,  k=%s, ωs=%s \n", rank, round(k, digits=6), round.(ωs[ik, :, 1]/physics.fz, digits=6))
    end    
end

counts = Nk_local * Nm * Neigs * (3 * grid.Ny + 1)
displs = cumsum(append!([0], counts))[1:comm_size]
MPI.Bcast!(counts, 0, comm)
MPI.Bcast!(displs, 0, comm)
MPI.Gatherv!(modes_local, rank == 0 ? VBuffer(modes, counts, displs) : nothing, 0, comm)

if rank == 0
    save_spectrum(ωs, modes, ks, ms, Neigs, grid.y, grid.Ny, files.nc)

    plot_growth_slice(files.nc, files.json, files.plotslicem)
    plot_growth(files.nc, files.json, files.plotgrowth)

    #plot_modes_1D(files.nc, files.json, files.plotmodes1D, Neigs)
    #plot_modes_2D(files.nc, files.json, files.plotmodes2D, Neigs)
end
