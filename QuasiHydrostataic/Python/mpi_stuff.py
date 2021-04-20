import numpy as np
from mpi4py import MPI             

def scatter_ks(Nk, dk, comm, rank, size):

    if rank == 0:
        # change Nk to Nk+1
        ks       = np.array(np.arange(0, Nk*dk, dk)).reshape(Nk,1)
        ks_split = np.array_split(ks, size, axis = 0)
        Nk_local = np.empty(0, int)

        for i in range(len(ks_split)):
            Nk_local = np.append(Nk_local, len(ks_split[i]))
            iks_ends = np.insert(np.cumsum(Nk_local), 0, 0)[0:-1]    
    else:
        ks       = None
        ks_split = None
        Nk_local = None
        iks_ends = None

    ks_split = comm.bcast(ks_split, root = 0)
    Nk_local = comm.bcast(Nk_local, root = 0)
    iks_ends = comm.bcast(iks_ends, root = 0)

    ks_local = np.zeros(np.shape(ks_split[rank]))
    comm.Scatterv([ks, Nk_local, iks_ends, MPI.DOUBLE], ks_local, root=0)

    return ks_local, Nk_local, iks_ends

def gather_ks(ks_local, Nk_local, Nk, comm, rank, size):

    if rank == 0:
        ks              = np.zeros((Nk,1), dtype=float)    
        Nk_local_input  = Nk_local
        Nk_local_output = Nk_local
        iks_ends_input  = np.insert(np.cumsum(Nk_local_input), 0,0)[0:-1]
        iks_ends_output = np.insert(np.cumsum(Nk_local_output),0,0)[0:-1]
    else:
        ks              = None
        Nk_local_input  = None
        Nk_local_output = None
        iks_ends_input  = None
        iks_ends_output = None

    Nk_local_output = comm.bcast(Nk_local_output, root = 0)
    iks_ends_output = comm.bcast(iks_ends_output, root = 0)

    comm.Barrier()

    comm.Gatherv(ks_local,[ks, Nk_local_output, iks_ends_output, MPI.DOUBLE], root=0)

    return ks

def gather_omegas(omegas_local, ks_local, Nk_local, Nk, Nm, Neigs, comm, rank, size):

    if rank == 0:
        omegas          = np.zeros((Nk, Nm, Neigs), dtype=complex)    
        Nk_local_input  = Nk_local * Nm * Neigs
        Nk_local_output = Nk_local * Nm * Neigs
        iks_ends_input  = np.insert(np.cumsum(Nk_local_input), 0,0)[0:-1]
        iks_ends_output = np.insert(np.cumsum(Nk_local_output),0,0)[0:-1]
    else:
        omegas          = None
        Nk_local_input  = None
        Nk_local_output = None
        iks_ends_input  = None
        iks_ends_output = None

    Nk_local_output = comm.bcast(Nk_local_output, root = 0)
    iks_ends_output = comm.bcast(iks_ends_output, root = 0)

    comm.Barrier()

    comm.Gatherv(omegas_local,[omegas, Nk_local_output, iks_ends_output, MPI.DOUBLE_COMPLEX], root=0)

    return omegas

def gather_modes(modes_local, ks_local, Nk_local, Nk, Nm, Neigs, Ny, comm, rank, size):

    if rank == 0:
        modes          = np.zeros((Nk, Nm, Neigs, 3*Ny+1), dtype=complex)    
        Nk_local_input  = Nk_local * Nm * Neigs * (3 * Ny + 1)
        Nk_local_output = Nk_local * Nm * Neigs * (3 * Ny + 1)
        iks_ends_input  = np.insert(np.cumsum(Nk_local_input), 0,0)[0:-1]
        iks_ends_output = np.insert(np.cumsum(Nk_local_output),0,0)[0:-1]
    else:
        modes           = None
        Nk_local_input  = None
        Nk_local_output = None
        iks_ends_input  = None
        iks_ends_output = None

    Nk_local_output = comm.bcast(Nk_local_output, root = 0)
    iks_ends_output = comm.bcast(iks_ends_output, root = 0)

    comm.Barrier()

    comm.Gatherv(modes_local,[modes, Nk_local_output, iks_ends_output, MPI.DOUBLE_COMPLEX], root=0)

    return modes
