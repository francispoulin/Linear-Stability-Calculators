# /usr/bin/env python

# LSA of QG jet using spectral (or finite differences)

# imports
import numpy as np
import scipy.linalg as spalg
import matplotlib.pyplot as plt
import sys


# Physical parameters
L     = 20.0                         # Length of domain
N     = 256                          # Number of grid points
F     = 0.1                          # 1/Bu = (L/Ld)^2
#beta  = 0.0                          # beta parameter

# Jet parameters
Lj = 1.0                           # width of jet
Uj = 1.0                           # maximum velocity of jet

method = 'cheb'


## FD2 computes the 2nd-order Finite Difference Operator
## ------
#    matrix on N+1 points (i.e. N intervals)
#    D2 = differentiation matrix
#    y  = uniform grid
def FD2(N,L):

    y  = np.linspace(0, L, N+1)       # Meridional coordiante
    dy = y[1] - y[0]                  # grid spacing
    dy2 = dy**2                       # square of dy

    Dyy = np.diag(1/dy2*np.ones(N),-1)   \
        + np.diag(-2/dy2*np.ones(N+1),0) \
        + np.diag(1/dy2*np.ones(N),1) 
    Dyy[0,0:3] = np.array([1, -2, 1]/dy2)
    Dyy[N,N-2:N+1] = np.array([1, -2, 1]/dy2)
    
    return y,Dyy


## CHEB computes the Chebyshev differentiation matrix
## ------
#    matrix on N+1 points (i.e. N intervals)
#    D = differentiation matrix
#    x = Chebyshev grid
def cheb(N):
    if N == 0:
        D = 0
        x = 1
    else:
        x = np.cos(np.pi*np.array(range(0,N+1))/N).reshape([N+1,1])
        c = np.ravel(np.vstack([2, np.ones([N-1,1]), 2])) \
            *(-1)**np.ravel(np.array(range(0,N+1)))
        c = c.reshape(c.shape[0],1)
        X = np.tile(x,(1,N+1))
        dX = X-(X.conj().transpose())
        D  = (c*(1/c).conj().transpose())/(dX+(np.eye(N+1)))   # off-diagonal entries
        D  = D - np.diag(np.sum(D,1))   # diagonal entries
    return D,x
## ------


# pick method
if method == 'FD2':
    ## FD2
    y,Dyy = FD2(N,L)
elif method == 'cheb':
    ## Cheb
    Dy,y = cheb(N)
    y   = (y[:,0]+1)*L/2
    Dy  = Dy*(2/L)
    Dyy = np.dot(Dy,Dy)
else:
    print("Method must be either FD2 or cheb")
    sys.exit()


# Define Basic State
P  = Uj*np.tanh((y-L/2)/Lj)
U  = Uj/(np.cosh((y-L/2)/Lj))**2
Q = - F*P + np.dot(Dyy,P) #+beta*y
Q2 = -2.*Uj*np.tanh((y-L/2)/Lj)/(np.cosh((y-L/2)/Lj))**2 
Qy = F*U - np.dot(Dyy,U) #+beta
I  = np.identity(N-1)

# Plot Velocity
plt.figure()
plt.clf()

plt.subplot(311)
plt.plot(y, P,'-b')
plt.grid(True)
plt.xlim([0, L])
plt.xlabel('y')
plt.title('Zonal Streamfunction')

plt.subplot(312)
plt.plot(y, U,'-b')
plt.grid(True)
plt.xlim([0, L])
plt.xlabel('y')
plt.title('Zonal Velocity')

plt.subplot(313)
plt.plot(y, Qy,'-b')
plt.grid(True)
plt.xlim([0, L])
plt.xlabel('y')
plt.title('Zonal PV_y')

plt.tight_layout()
plt.show()


# Define range of wavenumbers 
dk = 2e-2
kk = np.arange(dk,2+dk,dk)
Nk = len(kk)

# Define storage vectors
Ne = 2
grow = np.zeros((Ne,Nk))
freq = np.zeros((Ne,Nk))
modes = np.zeros((N+1,Ne,Nk),dtype=complex)

# Loop over wavenumbers
cnt = 0
for k in kk:

    # Useful Constants
    k2 = k**2

    # Set up Generalized Eigenvalue Problem (GEP)
    B = Dyy[1:-1,1:-1] - (k2 + F)*I
    A   = np.dot(np.diag(U[1:-1],0),B) + np.diag(Qy[1:-1],0) 
    
    # Solve for eigenvalues
    eigVals,eigVecs = spalg.eig(A,B)

    # Sort eigenvalues and eigenvectors
    ind = (-np.imag(eigVals)).argsort()
    eigVecs = eigVecs[:,ind]
    eigVals = k*eigVals[ind]

    # Store eigenvalues and eigenvectors
    grow[:,cnt] = eigVals[0:Ne].imag
    freq[:,cnt] = eigVals[0:Ne].real
    modes[1:N,:,cnt] = eigVecs[:,0:Ne]
    
    print ('Wavenumber (', int(cnt+1), '/', int(Nk),')', ': ',"{:.2f}".format(k),', Growth: ',"{:.4f}".format(grow[0,cnt]),', Phase: ',"{:.4f}".format(freq[0,cnt]))
    cnt += 1
    
print("Done!")


# Plot the two most unstable modes

plt.figure()
plt.clf()

# Initialize vector of maximum growth rate indices

Imax = np.zeros((2,1), dtype=int)

for ii in range(2):
    Imax[ii] = np.argmax(grow[ii,:])
    print("Max growth for curve", ii, "is", grow[ii,Imax[ii]])
    
plt.plot(kk,grow[0,:],'.b', label='Sinuous')
plt.plot(kk,grow[1,:],'.r', label='Varicos')
plt.grid(True)
plt.title("Growth Rates of QG (spectral) Jet: F = "+str(F))
plt.legend(loc='best')
plt.show()
#plt.savefig("QG_growth_bickley_spectral.png") 


# Plot Eigenvectors

plt.figure()
plt.clf()

for ii in range(2):
    plt.subplot(2,1,ii+1)
    plt.plot(y,modes[:,ii,Imax[ii]].real,'-b',linewidth=3,label=str("Mode "+str(ii+1)+": Re"))
    plt.plot(y,modes[:,ii,Imax[ii]].imag,'-r',linewidth=3,label=str("Mode "+str(ii+1)+": Im"))
    plt.legend()
    plt.grid(True)
    plt.title('Psi mode '+str(ii+1))
    
plt.tight_layout()
plt.show()


# Plot Mode Structure

plt.figure()
plt.clf()

for ii in range(2):
    Lx = 2*np.pi/kk[Imax[ii]]
    x = np.linspace(0,Lx,N+1)
    xx,yy = np.meshgrid(x,y)
    mode2d = (np.tile(modes[:,ii,Imax[ii]],(1,N+1))*np.exp(2*np.pi*1j*xx/Lx)).real
    plt.subplot(2,1,ii+1)
    plt.pcolormesh(xx/1e3,yy/1e3,mode2d,cmap='seismic')
    plt.xlim([0, Lx/1e3])
    plt.ylim([0, L/1e3])
    plt.colorbar()
    plt.xlabel('x (km)')
    plt.ylabel('y (km)')
    plt.title('Psi mode '+str(ii+1))
    
plt.tight_layout()
plt.show()    


