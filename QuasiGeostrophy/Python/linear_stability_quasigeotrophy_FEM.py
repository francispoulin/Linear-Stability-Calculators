# /usr/bin/env python

# LSA of QG jet using firedrake (FEM)

# Import Libraries

import numpy as np
import matplotlib.pyplot as plt

from firedrake import *
from firedrake.petsc import PETSc
try:
    from slepc4py import SLEPc
except ImportError:
    import sys
    warning("Unable to import SLEPc, eigenvalue computation not possible (try firedrake-update --slepc)")
    sys.exit(0)


# Define an Interval Mesh
Ly   = 20.0                     # domain width
Uj   = 1.0                      # jet strength
Lj   = 1.0                      # width of jet
n0   = 256                      # num grid pts
mesh = IntervalMesh(n0, Ly)     # compute grid 
x = SpatialCoordinate(mesh)     # spatial vect
yplot = np.linspace(0, Ly, n0)  # plotting vec
nPlotGrid = len(yplot)          # num plot pts

# Define parameters
#beta = Constant('0.1')          # f = f_0 + \beta y
F    = Constant('0.1')          # 1/Bu = (L/Ld)^2

# Profile
profile = 'bickley'
#profile = 'gaussian'
#profile = 'bump'

# order of the Method
p = 0 

# Define CG function space
V  = FunctionSpace(mesh,'CG',p+1)

# Impose zero Dirichlet BCs
bc = DirichletBC(V, 0.0, "on_boundary")
nodes = np.unique(bc.nodes)

# Translate grid as a vector
ydata = interpolate(x[0], V).vector().array()[1:]

# Define test and trial functions
phi, psi = TestFunction(V), TrialFunction(V)


# Background State

if profile == 'bickley':

    Pb  = Function(V).interpolate(Uj*tanh((x[0]-Ly/2)/Lj))
    Ub  = Function(V).interpolate(Uj/pow(cosh((x[0]-Ly/2)/Lj),2))
    dQb = Function(V).interpolate(4.0*Uj*(0.5 - pow(sinh((x[0]-Ly/2)/Lj),2))/pow(cosh((x[0]-Ly/2)/Lj),4) + F*(Uj/pow(cosh((x[0]-Ly/2)/Lj),2)))
    
elif profile == 'gaussian':
    
    Ub = Function(V).interpolate(Uj*exp(-(x[0]-Ly/2)**2))
    dQb = Function(V).interpolate(-2.*Uj*exp(-(x[0]-Ly/2)**2)*(2.*(x[0]-Ly/2)**2 - 1.))
    #dQb = Function(V).interpolate(-Uj*exp(-(x[0]-Ly/2)**2)*(4.*(x[0]-Ly/2)**2 - 4. - 1./(x[0]-Ly/2)**2))
    
elif profile == 'bump':
    
    Ub = Function(V).interpolate(conditional(lt(abs(x[0]-Ly/2.), 1.0), exp(-1./(1. -(x[0]-Ly/2)**2)), 0.0))
    dQb = Function(V).interpolate(conditional(lt(abs(x[0]-Ly/2.), 1.0), -2.*exp(-1./(1. -(x[0]-Ly/2)**2))*(3.*(x[0]-Ly/2.)**4 - 1)/((x[0]-Ly/2.)**2-1)**4, 0.0))

else:
    
    print("profile must be bickley, gaussian or bump")
    sys.exit()


# PLot Basic State

plt.figure()
plt.clf()

plt.subplot(311)
plt.plot(yplot,np.array([Pb.at(point, tolerance=1e-10) for point in yplot]),'-b')
plt.grid(True)
plt.xlabel('y')
plt.title('Zonal Streamfunction')

plt.subplot(312)
plt.plot(yplot,np.array([Ub.at(point, tolerance=1e-10) for point in yplot]),'-b')
plt.grid(True)
plt.xlabel('y')
plt.title('Zonal Velocity')

plt.subplot(313)
plt.plot(yplot,np.array([dQb.at(point, tolerance=1e-10) for point in yplot]),'-b')
plt.grid(True)
plt.xlabel('y')
plt.title('Zonal PV_y')

plt.tight_layout()
plt.show()


# set requested number of eigenvalues (use double the amount to account for signs. eg: want two modes, use four)
num_eigenvalues = 4

# Wavenumber
dk   = 2e-2
kk   = np.arange(dk, 2.+dk, dk)
kL   = len(kk)
egs_re  = np.zeros((len(kk),num_eigenvalues))
egs_im  = np.zeros((len(kk),num_eigenvalues))

# Initialize Eigenvector Storage
Eigenvalues = np.zeros((kL, 1, num_eigenvalues), dtype=complex)
Eigenvectors_re = np.zeros((kL, nPlotGrid, num_eigenvalues))
Eigenvectors_im = np.zeros((kL, nPlotGrid, num_eigenvalues))

# Define Petsc options
opts = PETSc.Options()
opts.setValue("eps_gen_non_hermitian", None)
opts.setValue("st_pc_factor_shift_type", "NONZERO")
opts.setValue("eps_type", "krylovschur")
opts.setValue("eps_largest_imaginary", None)
opts.setValue("eps_tol", 1e-10)

cnt = 0
for k in kk:

    # Define Useful Constants
    k2   = Constant(k**2)
    
    # Define Weak form
    a = ((Ub*phi).dx(0)*psi.dx(0))*dx + (((k2+F)*Ub - dQb)*phi*psi)*dx
    m = (phi.dx(0)*psi.dx(0) + (k2+F)*phi*psi)*dx
    petsc_a = assemble(a).M.handle
    petsc_m = assemble(m, bcs=bc).M.handle

    # Define Solver options
    es = SLEPc.EPS().create(comm=COMM_WORLD)
    es.setDimensions(num_eigenvalues)
    es.setOperators(petsc_a, petsc_m)
    es.setFromOptions()
    es.solve()

    # Additionally we can find the number of converged eigenvalues.
    nconv = es.getConverged()
    imax = min(nconv, num_eigenvalues)

    # Initialize List of Function objects (for debugging)
    #eval_real, eval_imag = [], []
    
    # Initialize Temporary Eigenvector Storage
    real, imag = Function(V), Function(V)
    
    # Find and Store Eigenvectors with their corresponding Eigenvalue 
    for i in range(imax):
        with real.dat.vec_wo as vr, imag.dat.vec_wo as vi:
            lam = es.getEigenpair(i, vr, vi)
        
        egs_re[cnt,i], egs_im[cnt,i] = k*lam.real, k*lam.imag
        
        #eval_real.append(real)
        #eval_imag.append(imag)
        #print("k = %f, i = %f, and eigenvalue = %f + 1j %f" % (k, i, lam.real, lam.imag))

        # Compute Eigenvectors as vectors to plot
        Eigenvectors_re[cnt,:,i] = np.array([real.at(point, tolerance=1e-12) for point in yplot])
        Eigenvectors_im[cnt, :,i] = np.array([imag.at(point, tolerance=1e-12) for point in yplot])

    print ('Wavenumber (', int(cnt+1), '/', int(kL),')', ': ',"{:.2f}".format(k),', Growth: ',"{:.4f}".format(egs_im[cnt,0]),', Phase: ',"{:.4f}".format(egs_re[cnt,0]))
    cnt += 1
    
print('Done!')
#print(np.max(abs(egs_im)))


for mode_index in [0,2]:
    ind = np.argmax(egs_im[:,mode_index])
    print('Mode:'+str(mode_index/2+1)+' Growth: '+str(egs_im[ind,mode_index]))

# Plot Growth Rates
plt.figure()
plt.plot(kk,egs_im[:,0],'.b',linewidth=2,label='sinuous')
plt.plot(kk,egs_im[:,2],'.r',linewidth=2,label='varicos')
plt.grid('on')
plt.xlabel('wavenumber')
plt.ylabel('growth')
plt.legend()
plt.title("Growth Rates of QG (firedrake) Jet: F = "+str(F.values()[0]))
plt.show()
#plt.savefig("QG_growth_bickley_firedrake.png") 


# Plot Eigenvectors
plt.figure()
plt.clf()

for mode_index in [0,2]:
    plt.subplot(2,1,int(mode_index/2+1))
    ind = np.argmax(egs_im[:,mode_index])
    plt.plot(yplot, Eigenvectors_re[ind,:,mode_index],'-b',linewidth=3, label='Re')
    plt.plot(yplot, Eigenvectors_im[ind,:,mode_index],'-r',linewidth=3, label='Im')
    plt.grid(True)
    plt.legend()
    plt.xlim([0, Ly])
    plt.title('Psi Mode:'+str(mode_index/2+1));
    
plt.tight_layout()    
plt.show()


# Plot Mode Structures
plt.figure()
plt.clf()

for mode_index in [0,2]:
    
    ind = np.argmax(egs_im[:,mode_index])
    psi_mode = Eigenvectors_re[ind,:,mode_index] + 1j*Eigenvectors_im[ind,:,mode_index]
    Lx = 2*np.pi/kk[ind]
    x = np.linspace(0,Lx,n0+1)
    xx,yy = np.meshgrid(x,ydata)
    psi_mode2d = (np.transpose(np.tile(psi_mode,(n0+1,1)))*np.exp(2*np.pi*1j*xx/Lx)).real
    
    plt.subplot(2,1,int(mode_index/2+1))
    plt.pcolormesh(xx,yy,psi_mode2d,cmap='seismic')
    plt.xlim([0, Lx])
    plt.ylim([0, Ly])
    plt.colorbar()
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Psi Mode '+str(mode_index/2+1));
    
plt.tight_layout()    
plt.show()



