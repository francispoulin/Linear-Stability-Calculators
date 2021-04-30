using MPI
using Printf

MPI.Init()

comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)
comm_size = MPI.Comm_size(comm)

# This assumes that N is divisible by comm_size!!!!
function find_split_size(test, comm_size, comm)

    N, M = size(test)

    counts = zeros(Int64, comm_size)
    displs = zeros(Int64, comm_size)

    counts[:] .= Int64(div(N, comm_size))
    displs[:]  = cumsum(append!([0], counts))[1:comm_size]

    return counts, displs
end

function array_split(test, counts, displs, comm_size, comm)

    N, M = size(test)

    split = zeros(Float64, comm_size, counts[1], M)
    for i in (1:comm_size)        
        split[i, :, :] = test[displs[i]+1:displs[i]+counts[i],:]
    end

    return split
end

N, M = 4, 3

if rank == 0
    test         = reshape(Float64.(1:M), 1, M)
    test         = repeat(test, N, 1)                              # from python example
    test         = permutedims(reshape((1.:(M*N)), M, N))          # a better test problem
    outputData_T = zeros((M, N))

    @printf("Original matrix: \n")
    @printf("================ \n")
    @printf("test is %s\n\n", test)

    counts, displs = find_split_size(test, comm_size, comm)      # need to generalize
 else
    test         = zeros(Float64, N, M)
    counts       = zeros(Int64, comm_size)
    displs       = zeros(Int64, comm_size)
    outputData_T = zeros((M, N))
end

test_T = permutedims(test)

MPI.Bcast!(counts, 0, comm)
MPI.Bcast!(displs, 0, comm)

if rank == 0
    split = array_split(test, counts, displs, comm_size, comm) 

    split_sizes_input  = counts*M
    split_sizes_output = counts*M

    displacements_input  = cumsum(append!([0], split_sizes_input))[1:comm_size]
    displacements_output = cumsum(append!([0], split_sizes_output))[1:comm_size]
    
    @printf("Scatter information:\n")
    @printf("====================\n")
    @printf("Input data split into vectors of sizes %s\n",       split_sizes_input)
    @printf("Input data split into displacements of sizes %s\n", displacements_input)
    @printf("\nSplit is of size %s\n\n", size(split))
else
    split = zeros(Float64, comm_size, counts[1], M)

    split_sizes_input    = zeros(Int64, comm_size)
    split_sizes_output   = zeros(Int64, comm_size)

    displacements_input  = zeros(Int64, comm_size)
    displacements_output = zeros(Int64, comm_size)
end

MPI.Bcast!(split, 0, comm)
MPI.Bcast!(split_sizes_output, 0, comm)
MPI.Bcast!(displacements_output, 0, comm)

output_chunk_T = permutedims(zeros(size(split[rank+1, :, :])))

MPI.Scatterv!(rank == 0 ? VBuffer(test_T, split_sizes_input, displacements_input) : nothing, output_chunk_T, 0, comm)

output = permutedims(output_chunk_T) #zeros(Float64, (size(output_chunk)[1], M))
if rank == 0
    @printf("Gathered array:\n")
    @printf("===============\n")
end
@printf("rank = %s  output = %s \n", rank, output)

MPI.Barrier(comm)

rank == 0 ? testnew = zeros(Float64, N) : nothing
MPI.Gatherv!(output_chunk_T, rank == 0 ? VBuffer(outputData_T, split_sizes_output, displacements_output) : nothing, 0, comm)
rank == 0 ? @printf("\nGathered outputData = %s\n", permutedims(outputData_T)) : nothing

MPI.Finalize()