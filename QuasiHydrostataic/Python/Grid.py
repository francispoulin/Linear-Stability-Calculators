import numpy as np

# FJP: must define Dy as well!!!!!

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

def Define_Grid(Ny, Ly, method):
    
    if   method == 'FD2':
        y,D2y = FD2(Ny, Ly)
    elif method == 'cheb':
        Dy,y = cheb(Ny)
        y   = y[:,0] * Ly / 2
        Dy  = Dy * ( 2 / Ly )
        D2y = np.dot(Dy,Dy)
    else:
        print("Method must be either FD2 or cheb")
        sys.exit()

    return y, Dy, D2y
