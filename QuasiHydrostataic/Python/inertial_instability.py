# /usr/bin/env python

"""
 Linear Stability Calculation of an inertially unstable jet using spectral methods

 The basic state is that of a barotropic Bickley jet (to be generalized).
 The perturbation is assumed to be periodic in both the zonal (x) and vertical (z) directions.

 Calculations are done in serial.  See inertial_instability_mpi.py for parallel execution.
"""
import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt
import sys

from src.Parameters import Grid, Physics, Jet, Files
from src.Parameters import Output_Parameters, save_spectrum
from src.Plotting_scripts import plot_growth_slice, plot_growth, plot_modes_1D, plot_modes_2D

### Define Parameters
file    = Files()
grid    = Grid(Ly = 1000e3, Lz = 3e3, Ny = 200, lat = np.pi/32)
physics = Physics(N=1e-2, nu=0.26, kwargs={"lat": grid.lat, "NT": 1})
jet     = Jet(kwargs={"y": grid.y, "Ly": grid.Ly, "fz": physics.fz, "fy": physics.fy})

### Output Parameters 
Output_Parameters(grid, physics, jet, file.json)

### Number of eigenvalues to store and wavenumbers
Neigs  = 10                            
dk, Nk = 1e-6, 40
dm, Nm = 1e-4, 150

iks = np.array(range(Nk)) 
ims = np.array(range(1,Nm+1)) 
ks  = np.array([dk*ik for ik in iks])      
ms  = np.array([dm*im for im in ims])      

omegas = np.zeros((Nk, Nm, Neigs),              dtype=complex)
modes  = np.zeros((Nk, Nm, Neigs, 3*grid.Ny+1), dtype=complex)

Z, I     = np.zeros((grid.Ny+1, grid.Ny+1)), np.eye(grid.Ny+1)
Qm, Bym  = np.diag(physics.fz - jet.dU), np.diag(- physics.fy * jet.dU)
N2       = physics.N**2

### Loop over wavenumbers
for (ik, k) in enumerate(ks):
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
            omegas[ik, im, ie] = eigVals[ie]
            modes [ik, im, ie] = eigVecs[:,ie]

save_spectrum(omegas, modes, ks, ms, jet.y, Neigs, grid.Ny, file.nc)

plot_growth_slice(file.nc, file.json, file.plotslicem)
plot_growth(file.nc, file.json, file.plotgrowth)
plot_modes_1D(file.nc, file.json, file.plotmodes1D, Neigs)
plot_modes_2D(file.nc, file.json, file.plotmodes2D, Neigs)

# To-Do
# -> create: build_A, compute_spectrum, sort perturbation 
# -> better organize files
# -> Julia
# -> 2D eigenvalue problem
