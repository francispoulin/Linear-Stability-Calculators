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


function clencurt(N)

     θ = π * collect(0:N) ./N
     x = cos.(θ)
     w = zeros(1, N+1)
    ii = collect(2:N)
     v = ones(N-1,1)
    
    if mod(N,2) == 0
          w[1] = 1/(N^2 - 1)
        w[N+1] = w[1]
        for k=1:N/2-1
             v = v .- 2 * cos.(2 * k * θ[ii])/(4 * k^2 - 1)
        end
        v = v - cos.(N*θ[ii]) / (N^2 - 1)
    else
          w[1] = 1/N^2
        w[N+1] = w[1]
        for k=1:(N-1)/2
            v = v .- 2 * cos.(2 * k * θ[ii])/(4 * k^2 - 1)
        end
    end
    
    w[ii] = 2*v/N

    return x, w
end
