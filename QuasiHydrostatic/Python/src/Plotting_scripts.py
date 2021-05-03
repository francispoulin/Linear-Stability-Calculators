import numpy as np
import scipy.linalg as spalg
import matplotlib.pyplot as plt

def plot_profiles(L, y, U, Q, B, By, file):
    
    plt.clf()

    fig, axs = plt.subplots(4,1)
    fig.suptitle('Basic State Profiles')

    axs[0].plot(y / L, U,  '-b',linewidth=2, label='U')
    axs[1].plot(y / L, Q,  '-b',linewidth=2, label='Q')
    axs[2].plot(y / L, B,  '-b',linewidth=2, label='B') 
    axs[3].plot(y / L, By, '-b',linewidth=2, label='By')

    for (i, ax) in enumerate(fig.get_axes()):
        ax.label_outer()
        ax.legend()
        axs[i].grid('on')

    axs[3].set(xlabel='y (km)')
    plt.savefig(file)

import xarray as xr
import json

def plot_growth_slice(filenc, filejson, fileplot):

    # read json file
    with open(filejson) as f:
        data = json.load(f)

    Umax = data["jet"]["Umax"]
    L    = data["jet"]["L"]
    fz   = data["physics"]["fz"]

    # read NetCDF file
    ds = xr.open_dataset(filenc)

    growths = ds["omegas_imag"]
    Neigs = len(growths[0,0,:])
    ms = ds["m"]

    print("\n--> Plotting the growth rates for a fixed k versus m in ", fileplot)

    # make plot
    plt.clf()
    for i in range(Neigs):
        plt.plot(L/(2*np.pi)*ms, growths[0, :, i]/fz) 
        
        plt.xlabel('mL/2*π')
        plt.xlim([L/(2*np.pi)*ms[0], L/(2*np.pi)*ms[-1]])
        plt.ylim([0, 1.1*np.amax(growths)/fz])
        plt.grid('on')
        plt.title('Growth Rates')
        plt.savefig(fileplot)

def plot_growth(filenc, filejson, fileplot):

    # read json file
    with open(filejson) as f:
        data = json.load(f)

    Umax = data["jet"]["Umax"]
    L    = data["jet"]["L"]
    fz   = data["physics"]["fz"]

    # read NetCDF file
    ds = xr.open_dataset(filenc)

    growths = ds["omegas_imag"]
    Neigs = len(growths[0,0,:])
    ks = ds["k"]
    ms = ds["m"]

    print("\n--> Plotting the growth rates versus (k, m) in ", fileplot)

    # make plot
    plt.clf()
    plt.contourf(L/(2*np.pi)*ks, 
                 L/(2*np.pi)*ms, 
                (np.maximum(0,growths[:,:,0])/fz).T)
    plt.xlabel('kL/2*π')
    plt.ylabel('mL/2*π')
    plt.colorbar()
    plt.title('First most unstable mode')
    plt.savefig(fileplot)

def plot_modes_1D(filenc, filejson, plotmodes1D, Neigs):

    # read json file
    with open(filejson) as f:
        data = json.load(f)

    L      = data["jet"]["L"]
    Ny     = data["grid"]["Ny"]
    Ly     = data["grid"]["Ly"]
    method = data["grid"]["method"]
    fz     = data["physics"]["fz"]  

    # read NetCDF file
    ds = xr.open_dataset(filenc)

    modes_real = ds["modes_real"]
    modes_imag = ds["modes_imag"]
    growth     = ds["omegas_imag"]
    k          = ds["k"]
    m          = ds["m"]
    y          = ds["ys"][0:Ny+1]

    ik   = 0                   # since this is where we find the most unstable mode

    print("\n--> Plotting 1D modes in files with the name", plotmodes1D)
    for iEig in range(Neigs):
        im   = growth[ik, :, iEig].argmax()

        print("     iEig = ", iEig,
              "k = % 10.6f" %    float(L/(2*np.pi)*k[ik]),
              "m = % 10.6f" %      float(L/(2*np.pi)*m[im]),
              "growth = % 10.6f" % float((growth[ik, im, iEig]/fz)))

        vvec = np.zeros(Ny+1, dtype=complex)
        uvec       = modes_real[ik, im, iEig, 0:Ny+1]      + 1j * modes_imag[ik, im, iEig, 0:Ny+1]
        vvec[1:-1] = modes_real[ik, im, iEig, Ny+1:2*Ny]   + 1j * modes_imag[ik, im, iEig, Ny+1:2*Ny]
        bvec       = modes_real[ik, im, iEig, 2*Ny:3*Ny+1] + 1j * modes_imag[ik, im, iEig, 2*Ny:3*Ny+1]

        plt.clf()
        fig, axs = plt.subplots(3,1)
        fig.suptitle('Perturbation Profiles')

        axs[0].plot(y / L, uvec.real,  '-b',linewidth=2, label='Re(u)')
        axs[0].plot(y / L, uvec.imag,  '-r',linewidth=2, label='Im(u)')
        axs[1].plot(y / L, vvec.real,  '-b',linewidth=2, label='Re(v)')
        axs[1].plot(y / L, vvec.imag,  '-r',linewidth=2, label='Im(v)')
        axs[2].plot(y / L, bvec.real,  '-b',linewidth=2, label='Re(b)') 
        axs[2].plot(y / L, bvec.imag,  '-r',linewidth=2, label='Re(b)') 

        for (i, ax) in enumerate(fig.get_axes()):
            ax.label_outer()
            ax.legend()
            axs[i].grid('on')

        axs[2].set(xlabel='y (km)')
        file = 'modes_1D_iEig' + str(iEig) + '.png'
        plt.savefig(file)
        plt.close()


def plot_modes_2D(filenc, filejson, plotmodes2D, Neigs):

    # read json file
    with open(filejson) as f:
        data = json.load(f)

    L      = data["jet"]["L"]
    Ny     = data["grid"]["Ny"]
    Ly     = data["grid"]["Ly"]
    method = data["grid"]["method"]
    fz = data["physics"]["fz"]  

    # read NetCDF file
    ds = xr.open_dataset(filenc)

    modes_real = ds["modes_real"]
    modes_imag = ds["modes_imag"]
    growth     = ds["omegas_imag"]
    k          = ds["k"]
    m          = ds["m"]
    y          = ds["ys"][0:Ny+1]

    ik   = 0
    iEig  = 0
    im    = growth[ik, :, iEig].argmax()
    Nz    = Ny+1
    Lz    = 2*np.pi/m[im]
    z     = np.arange(0, Lz, Lz/Nz)
    [Y,Z] = np.meshgrid(y,z)

    print("\n--> Plotting 2D modes in files with the name", plotmodes2D)
    for iEig in range(Neigs):

        im   = growth[ik, :, iEig].argmax()

        print("     iEig = ", iEig,
              "k = % 10.6f" %    float(L/(2*np.pi)*k[ik]),
              "m = % 10.6f" %      float(L/(2*np.pi)*m[im]),
              "growth = % 10.6f" % float((growth[ik, im, iEig]/fz)))

        vvec       = np.zeros(Ny+1, dtype=complex)
        uvec       = modes_real[ik, im, iEig, 0:Ny+1]      + 1j * modes_imag[ik, im, iEig, 0:Ny+1]
        vvec[1:-1] = modes_real[ik, im, iEig, Ny+1:2*Ny]   + 1j * modes_imag[ik, im, iEig, Ny+1:2*Ny]
        bvec       = modes_real[ik, im, iEig, 2*Ny:3*Ny+1] + 1j * modes_imag[ik, im, iEig, 2*Ny:3*Ny+1]

        u = np.tile(uvec.real,(Nz,1)) * np.cos(float(m[im])*Z) - np.tile(uvec.imag,(Nz,1)) * np.sin(float(m[im])*Z)
        v = np.tile(vvec.real,(Nz,1)) * np.cos(float(m[im])*Z) - np.tile(vvec.imag,(Nz,1)) * np.sin(float(m[im])*Z)
        b = np.tile(bvec.real,(Nz,1)) * np.cos(float(m[im])*Z) - np.tile(bvec.imag,(Nz,1)) * np.sin(float(m[im])*Z)

        plt.clf()
        fig, axs = plt.subplots(1, 3, figsize=(30,10))
        fig.suptitle('Perturbation Profiles: 2D', fontsize=20)

        axs[0].pcolormesh(Y/1e3, Z/1e3, u, shading='gouraud')
        axs[1].pcolormesh(Y/1e3, Z/1e3, v, shading='gouraud')
        axs[2].pcolormesh(Y/1e3, Z/1e3, b, shading='gouraud')
        axs[0].set_title('u', fontsize=20)
        axs[1].set_title('v', fontsize=20)
        axs[2].set_title('b', fontsize=20)

        for (i, ax) in enumerate(fig.get_axes()):
            ax.label_outer()
            ax.set_xlabel('y (km)')

        axs[0].set(ylabel='z (km)')
        file = 'modes_2D_iEig' + str(iEig) + '.png'
        plt.savefig(file)
        plt.close()
