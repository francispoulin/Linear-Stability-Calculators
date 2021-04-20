import numpy as np
from dataclasses import dataclass, asdict, field
from typing import Any

from Grid import Define_Grid
from Plotting_scripts import plot_profiles

from dataclasses import asdict
import json

def without_keys(d, ex): 
    for i in ex: d.pop(i) 
    return d

@dataclass
class Files:
    """Class to store the different file names. """
    nc:          str = 'data_spectrum.nc'
    json:        str = 'data_parameters.json'
    plotgrowth:  str = 'growth.png'
    plotslicem:  str = 'growth_slice_m.png'
    plotmodes1D: str = 'growth_modes1D'
    plotmodes2D: str = 'growth_modes2D'

@dataclass
class Grid:
    """Class for storing grid parameters."""
    Ly:      float = 1000e3
    Lz:      float = 3e3
    Ny:      int   = 100
    lat:     float = np.pi/4
    method:  str   = 'cheb'

    y:       Any   = field(init=False)
    Dy:      Any   = field(init=False)
    D2y:     Any   = field(init=False)

    def __post_init__(self):
        self.y, self.Dy, self.D2y = Define_Grid(self.Ny, self.Ly, self.method)

@dataclass
class Physics:
    """Class for storing physical parameters. """
    Omega:   float = 2*np.pi/(24*3600)
    N:       float = 1e-3
    g:       float = 9.81
    nu:      float = 0.26

    fz:      float = field(init=False)        # change to any if beta-plane
    fy:      float = field(init=False)        # change to any if beta-plane

    kwargs: field(default_factory=dict) = None

    def __post_init__(self):
        [setattr(self, k, v) for k, v in self.kwargs.items()]

        self.fz = 2*self.Omega * np.sin( self.lat )
        self.fy = 2*self.Omega * np.cos( self.lat ) * self.NT

@dataclass
class Jet:
    """Class for the basic state profile. """
    L:        float = field(init=False) 
    Umax:     float = 14.6 
    profile:  str   = 'Bickley'

    iU:        Any   = field(init=False)
    U:         Any   = field(init=False)
    dU:        Any   = field(init=False)

    kwargs: field(default_factory=dict) = None

    def __post_init__(self):      
        [setattr(self, k, v) for k, v in self.kwargs.items()]

        self.L = self.Ly/10
        self.iU, self.U, self.dU = Set_Profile(self.y, self.L, self.Umax, self.fz, self.fy, self.profile)

def Output_Parameters(grid, physics, jet, file):

    fz = physics.fz
    N  = physics.N
    g  = physics.g
    nu = physics.nu

    U = jet.Umax
    L = jet.L

    # Non-dimensional parameters
    phi0  = U * fz * L
    Ro    = U / (fz * L)
    Bu    = (N / fz)**2
    delt  = phi0/(g * L)
    gamma = g**2 * nu / (fz*phi0**2)

    # FJP: Is this the output you really want???
    # Output parameters in terminal
    print(' ')
    print('Nondimensional Parameters')
    print('=========================')
    print('Ro      = ', Ro)
    print('Bu      = ', Bu)
    print('phi0    = ', phi0)
    print('delta   = ', delt)
    print('gamma   = ', gamma)
    print('method  = ', grid.method)
    print('profile = ', jet.profile)
    print(' ')

    # Output parameters to a json file
    json_data = {
    'grid':    without_keys(asdict(grid),['y', 'Dy', 'D2y']),
    'physics': without_keys(asdict(physics),['kwargs']),
    'jet':     without_keys(asdict(jet),['iU', 'U', 'dU', 'kwargs'])
    }

    with open(file, 'w') as fp:
        json.dump(json_data, fp, sort_keys=True, indent=4)

import xarray as xr 

def save_spectrum(omegas, modes, ks, ms, y, Neigs, Ny, file):

    print("\n--> Saving the spectrum and modal structures in ", file)

    ys = np.hstack([y, y[1:-1], y])
    
    ds = xr.Dataset(
        {
            "omegas_real": (("k","m","i"),      omegas.real),
            "omegas_imag": (("k","m","i"),      omegas.imag),
            "modes_real":  (("k","m","i","ys"), modes.real),
            "modes_imag":  (("k","m","i","ys"), modes.imag),
            },
        coords={
            "k":  ks.flatten(),
            "m":  ms,
            "i":  range(Neigs),
            "ys": ys,
        },
    )

    ds.to_netcdf(file)


def Set_Profile(y, L, Umax, fz, fy, profile):

    if profile == 'Bickley':
        iU, U, dU = Bickley(y, L, Umax, fz, fy)
    else:
        print("profile must be Bickley")
        sys.exit()

    return iU, U, dU

"""
Compute the velocity profile to be a Bickley jet
"""
def Bickley(y, L, Umax, fz, fy):

    iU  =       L * Umax *  np.tanh( y / L )
    U   =           Umax / (np.cosh( y / L )**2 )
    dU  =   2 / L * Umax *  np.tanh( y / L )/(np.cosh( y / L )**2)
    Q   =   fz - dU
    B   = - fy * U
    By  = - fy * dU

    plot_profiles(L, y, U, Q, B, By, 'structure_jet.png')

    return iU, U, dU
