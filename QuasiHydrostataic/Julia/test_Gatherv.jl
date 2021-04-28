using MPI
using Printf

MPI.Init()

comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)
size = MPI.Comm_size(comm)

function array_split(N, size, comm)

    counts = zeros(Int64, size)
    displs = zeros(Int64, size)

    if rank == 0
        counts_guess = Int64(div(N, size, RoundDown))
        Remainder    = Int64(N - counts_guess*size)
        counts[:]   .= counts_guess

        for i in 1:Remainder
            counts[i] += 1
        end

        displs[:] = cumsum(append!([0], counts))[1:size]
    end

    MPI.Bcast!(counts, 0, comm)
    MPI.Bcast!(displs, 0, comm)

    return counts, displs
end

N = 10

### Find counts and displacement arrays for scattering
counts, displs = array_split(N, size, comm)

### Define array
rank == 0 ? k = Float64[i for i in 1:N] : k = nothing
rank == 0 ? @printf("original k = %s \n", k) : nothing

### Scatter array 
k_local = zeros(Float64, counts[rank+1])
MPI.Scatterv!(rank == 0 ? VBuffer(k, counts, displs) : nothing, k_local, 0, comm)
@printf("rank = %s  k_local = %s \n", rank, k_local)

### Gatter array
rank == 0 ? knew = zeros(Float64, N) : nothing
MPI.Gatherv!(k_local, rank == 0 ? VBuffer(knew, counts, displs) : nothing, 0, comm)
rank == 0 ? @printf("gathered k = %s\n", knew) : nothing

MPI.Finalize()
