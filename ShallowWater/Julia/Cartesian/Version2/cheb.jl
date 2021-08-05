using LinearAlgebra

function cheb(N)
    if N==0 
        D=0; x=1
        return D,x
    else
        x= cos.(pi * collect(0:N) ./ N)
        c=[2; ones(N-1, 1); 2] .* (-1).^collect(0:N);
        X = repeat(x, 1, N+1);
        dX=X .- X';
        D = (c .* (1 ./ c)')  ./ (dX .+ Diagonal( ones(N+1) ) );
        D = D .- Diagonal( vec( sum(D', dims=1) ) );
        return D,x
    end
end