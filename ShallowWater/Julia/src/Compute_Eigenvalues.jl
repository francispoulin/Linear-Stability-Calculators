using LinearAlgebra

"""
Cartesian version
"""
function build_matrix(geometry::Cartesian;
    k,
    solution, 
    phys
    )

    N = phys.grid.N
    g = phys.g

    k2 = k^2
    ik = 1im*k

    u = solution.u
    h = solution.η .+ phys.grid.H
    D = phys.grid.D

     U = diagm(0 => u)
     H = diagm(0 => h)
    dU = diagm(0 => D*u)
    dH = diagm(0 => D*h)
     F = diagm(0 => phys.coriolis)
     I = diagm(0 => ones(N+1))

    # Form 1L-RSW Matrix
    #   [u1, v1, h1]

    A = [  U            (-F+dU)[:, 2:N]      g*I;
          -F[2:N,:]/k2       U[2:N,2:N]  -g/k2*D[2:N,:];
           H               D*H[:,  2:N]        U];

    return A 
end

"""
Spherical version
"""
function build_matrix(geometry::Spherical;
    k,
    solution,
    phys
    )

    a = phys.grid.a
    N = phys.grid.N
    g = phys.g
    
    k2 = k^2
    ik = 1im*k    

        u = solution.u
        h = solution.η .+ phys.grid.H
        D = phys.grid.D
        ϕ = phys.grid.ϕ
    axcos = a * cos.(ϕ)
    
           D = D/a
       IxCos = diagm(0=>axcos)
       IoCos = diagm(0=>1 ./ axcos)
    TwoΩxSin = diagm(0=>phys.coriolis);    #
       UoCos = diagm(0=>u./axcos)
       UxTan = diagm(0=>u./axcos.*sin.(ϕ))
       HoCos = diagm(0=>h./axcos)
       HxCos = diagm(0=>h.*axcos)
          dU = diagm(0=>D*u);
          dH = diagm(0=>D*h);

    # Form 1L-RSW Matrix
    #   [u1, v1, h1]

    A = [                UoCos             (-TwoΩxSin - UxTan + dU)[:, 2:N]    g*IoCos;
          -(TwoΩxSin + 2*UxTan)[2:N,:]/k2                     UoCos[2:N,2:N]  -g/k2*D[2:N,:];
                         HoCos                  (IxCos * D * HxCos)[:, 2:N]     UoCos];

    return A
end

function compute_spectrum(A, k, phys, Nmodes)

    λ_in, λvec_in = eigen(A);
            index = sortperm(imag(λ_in), rev=true);
             λvec = λvec_in[:,index];
                λ = λ_in[index];
     
           growth = imag(λ[1:Nmodes]) * k;
             freq = real(λ[1:Nmodes]) * k;
            modes = λvec[:,1:Nmodes];
    
        return growth, freq, modes
    end