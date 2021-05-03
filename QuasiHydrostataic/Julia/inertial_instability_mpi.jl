#!/usr/bin/env julia --startup-file=no

"""
Linear-Stability-Calculators: Quasi-Hydrostatic Model
=====================================================

Linear Stability Calculation of an inertially unstable jet using spectral methods

The basic state is that of a barotropic Bickley jet (to be generalized).
The perturbation is assumed to be periodic in both the zonal (x) and vertical (z) directions.

Calculations are done in parallel.  See inertial_instability.py for serial execution.
"""

# put this somewhere else?
"""
    split_count(N::Integer, n::Integer)
Return a vector of `n` integers which are approximately equally sized and sum to `N`.
"""
function split_count(N::Integer, n::Integer)
    q,r = divrem(N, n)
    return [i <= r ? q+1 : q for i = 1:n]
end

using MPI
using Printf

MPI.Init()

comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)
comm_size = MPI.Comm_size(comm)

# put this elsewhere?
abstract type AbstractMesh end

struct Cheb <: AbstractMesh end
struct FiniteDifferenceOrder2  <: AbstractMesh end

include("src/parameters.jl")
include("src/mpi_stuff.jl")
include("src/Compute_Eigenvalues.jl")
include("src/Plotting_Functions.jl")

files   = Files()
grid    = Grid(Ly = 1000e3, Lz = 3e3, Ny = 300, θ₀ = π/32, method=Cheb())
physics = Physics(N = 1e-2, ν = 0.26, θ₀ = grid.θ₀, NT = 1)
jet     = Jet(grid, physics)

rank == 0 ? Output_Parameters(grid, physics, jet, files.json) : nothing

Neigs  = 10                            
dk, Nk = 1e-6, 41
dm, Nm = 1e-4, 150

ms = collect(dm:dm:Nm*dm)

#rank == 0 ? ks = Float64[(ik-1)*dk for ik in 1:Nk] : ks = nothing

if rank == 0
    ks = Float64[(ik-1)*dk for ik in 1:Nk]
    Nk_counts = split_count(Nk, comm_size)
    ks_vbuf   = VBuffer(ks, Nk_counts)        # VBuffer for scatter

else
    Nk_counts = nothing
    ks_vbuf   = nothing
end

Nk_local, = MPI.Scatter!(Nk_counts, zeros(Int, 1), 0, comm)
ks_local = MPI.Scatterv!(ks_vbuf, zeros(Float64, Nk_local), 0, comm)

#@printf("rank = %s, Nk_local = %s, ks_local = %s\n", rank, Nk_local, ks_local)

ωs_local    = zeros(Complex{Float64}, (Neigs, Nm, length(ks_local)))
modes_local = zeros(Complex{Float64}, (3*grid.Ny+1, Neigs, Nm, length(ks_local)))

for (ik, k) in enumerate(ks_local)
    for (im, m) in enumerate(ms)
        A = build_matrix(k, m, grid, physics, jet)

        ωs_local[:, im, ik], modes_local[:, :, im, ik] = compute_spectrum(A, Neigs)

        print("k = ", ks_local[ik], " m = ", ms[im], " growth = ", imag(ωs_local[1, im, ik]), "\n")
    end
end

MPI.Barrier(comm)

if rank == 0
    ωs    = zeros(Complex{Float64}, (Neigs, Nm, Nk)) 
    modes = zeros(Complex{Float64}, (3*grid.Ny+1, Neigs, Nm, Nk))

    Neigs_counts = [Neigs for k = 1:comm_size]
    Nm_counts    = [Nm    for j = 1:comm_size]

    ωs_sizes  = vcat(Neigs_counts', Nm_counts', Nk_counts')
    ωs_counts = vec(prod(ωs_sizes, dims=1))
    ωs_vbuf   = VBuffer(ωs, ωs_counts)   

    Ny_counts = [3*grid.Ny+1 for k = 1:comm_size]

    modes_sizes  = vcat(Ny_counts', Neigs_counts', Nm_counts', Nk_counts')
    modes_counts = vec(prod(modes_sizes, dims=1))
    modes_vbuf   = VBuffer(modes, modes_counts) 

else
    ωs_sizes = nothing
    ωs_vbuf  = nothing    

    modes_sizes = nothing
    modes_vbuf  = nothing    
end

MPI.Gatherv!(ωs_local,    ωs_vbuf,    0, comm)
MPI.Gatherv!(modes_local, modes_vbuf, 0, comm)

if rank == 0
    @printf("ωs size = %s, modes size = %s \n", size(ωs), size(modes))
end
#if rank == 0
#    for (ik, k) in enumerate(ks)
#        @printf("r=%s,  k=%s, ωs=%s \n", rank, round(k, digits=6), round.(ωs[ik, :, 1]/physics.fz, digits=6))
#    end    
#end

if rank == 0
    save_spectrum(ωs, modes, ks, ms, Neigs, grid.y, grid.Ny, files.nc)

    plot_growth_slice(files.nc, files.json, files.plotslicem)
    plot_growth(files.nc, files.json, files.plotgrowth)

    plot_modes_1D(files.nc, files.json, files.plotmodes1D, Neigs)
    plot_modes_2D(files.nc, files.json, files.plotmodes2D, Neigs)
end

MPI.Finalize()
