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
