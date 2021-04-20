# /usr/bin/env python

"""
 Linear Stability Calculation of an inertially unstable jet using spectral methods

 The basic state is that of a barotropic Bickley jet (to be generalized).
 The perturbation is assumed to be periodic in both the zonal (x) and vertical (z) directions.

"""
import numpy as np
from numpy import linalg as LA
from mpi4py import MPI             
import matplotlib.pyplot as plt
import sys

from Parameters import Grid, Physics, Jet, Files
from Parameters import Output_Parameters, save_spectrum

from Plotting_scripts import plot_growth_slice, plot_growth, plot_modes_1D, plot_modes_2D

from mpi_stuff import scatter_ks, gather_spectrum

comm = MPI.COMM_WORLD              
rank = comm.Get_rank()             
size = comm.Get_size()             

### Define Parameters
file    = Files()
grid    = Grid(Ly = 1000e3, Lz = 3e3, Ny = 10, lat = np.pi/32)
physics = Physics(N=1e-2, nu=0.26, kwargs={"lat": grid.lat, "NT": 1})
jet     = Jet(kwargs={"y": grid.y, "Ly": grid.Ly, "fz": physics.fz, "fy": physics.fy})

### Output Parameters 
if rank==0:
    Output_Parameters(grid, physics, jet, file.json)

### Number of eigenvalues to store and wavenumbers
Neigs  = 10                            
dk, Nk = 1e-6, 10 #40
dm, Nm = 1e-4, 20 #150
ms     = np.arange(dm, Nm*dm, dm)

ks_local, Nk_local, iks_ends = scatter_ks(Nk, dk, comm, rank, size)

omegas_local = np.zeros((Nk_local[rank], Nm, Neigs),              dtype=complex)
modes_local  = np.zeros((Nk_local[rank], Nm, Neigs, 3*grid.Ny+1), dtype=complex)

Z, I     = np.zeros((grid.Ny+1, grid.Ny+1)), np.eye(grid.Ny+1)
Qm, Bym  = np.diag(physics.fz - jet.dU), np.diag(- physics.fy * jet.dU)
N2       = physics.N**2

### Loop over wavenumbers
for (ik, k) in enumerate(ks_local):
    k2 = k**2
    for (im, m) in enumerate(ms):
        m2 = m**2

        Um    = k*np.diag(jet.U,0) - 1j*physics.nu*m2*I
        fyomDy = physics.fy/m*grid.Dy

        # Set up eigenvalue problem: [u, v, b]
        A = np.vstack((np.hstack((                              Um,                  (1j*Qm + fyomDy)[:,1:-1],     - 1j*k/m*I)),
                       np.hstack(( - (1j * physics.fz * I + fyomDy)[1:-1,:],                    Um[1:-1,1:-1],   - 1/m*grid.Dy[1:-1,:])),
                       np.hstack((                     1j*N2*k/m*I,          (N2/m*grid.Dy - 1j* Bym)[:,1:-1],             Um)) ))

        # Compute eigenspectrum
        eigVals,eigVecs = LA.eig(A)

        # Sort eigenvalues
        ind = (-np.imag(eigVals)).argsort()
        eigVecs = eigVecs[:,ind]
        eigVals = eigVals[ind]

        # Store spectrum
        for ie in range(Neigs):
            if rank == 0:
                omegas_local[ik, im, ie] = eigVals[ie]
                modes_local [ik, im, ie] = eigVecs[:,ie]

omegas_real, omegas_imag = gather_spectrum(omegas_local, ks_local, Nk_local, Nk, Nm, Neigs, comm, rank, size)

print("rank = ", rank, " and shape = ", omegas_real.shape)
sys.exit()

# To-Do
# -> assign complex omegas array correctly
# -> generalize gather_spectrum to work for ks (real), omegas (complex), modes (complex)

omegas = np.array((Nk, Nm, Neigs), dtype=complex)
omegas = omegas_real #+ 1j*omegas_imag

# Gather to 0
ks     = np.array(comm.gather(ks_local,     root=0))   
#omegas = np.array(comm.gather(omegas, root=0))   
modes  = np.array(comm.gather(modes_local,  root=0))   

if rank == 0:                                    
    ks     = np.reshape(ks,     Nk)                      
#    omegas = np.reshape(omegas, (Nk, Nm, Neigs))                      
    modes  = np.reshape(modes,  (Nk, Nm, Neigs, 3*grid.Ny+1))                      

if rank == 0:
    save_spectrum(omegas, modes, ks, ms, jet.y, Neigs, grid.Ny, file.nc)

    plot_growth_slice(file.nc, file.json, file.plotslicem)
    plot_growth(      file.nc, file.json, file.plotgrowth)
    plot_modes_1D(    file.nc, file.json, file.plotmodes1D, Neigs)
    plot_modes_2D(    file.nc, file.json, file.plotmodes2D, Neigs)

# To-Do
# -> numba parallel k loop using krange!
# -> create src directory 
# -> create: build_A, compute_spectrum, sort perturbation 
# -> better organize files
# -> Julia
# -> 2D eigenvalue problem
