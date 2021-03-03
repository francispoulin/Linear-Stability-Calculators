using LinearAlgebra

"""
Cartesian version of build_matrix
"""
function build_matrix(geometry::Cartesian;
    k,
    background, 
    phys
    )

    solution = background[geometry]

    N = phys[geometry].grid.N
    g = phys[geometry].g

    kx2 = k^2
    ik = 1im*k

     u = solution.u
     h = solution.η .+ phys[geometry].grid.H
    Dy = phys[geometry].grid.D

     U   = diagm(0 => u)
     H   = diagm(0 => h)
    dUdy = diagm(0 => Dy*u)
     F   = diagm(0 => phys[geometry].coriolis)
     I   = diagm(0 => ones(N+1))

    # Form 1L-RSW Matrix
    #   [u1, v1, h1]

    A = [  U             (-F+dUdy)[:, 2:N]         g*I;
          -F[2:N,:]/kx2          U[2:N,2:N]  -g/kx2*Dy[2:N,:];
           H                  Dy*H[:,  2:N]          U];

    return A 
end

"""
Spherical version of build_matrix
"""
function build_matrix(geometry::Spherical;
    k,
    background,
    phys
    )

    solution = background[geometry]

    a = phys[geometry].grid.a
    N = phys[geometry].grid.N
    g = phys[geometry].g
    
    kλ2 = k^2
    ik = 1im*k    

         u = solution.u
         h = solution.η .+ phys[geometry].grid.H
        Dy = phys[geometry].grid.D / a          
         ϕ = phys[geometry].grid.ϕ
    axcos = a * cos.(ϕ)
    
       IoACos = diagm(0=>1 ./ axcos)
            F = diagm(0=>phys[geometry].coriolis);    #
       UoACos = diagm(0=>u./axcos)
      UxTanoA = diagm(0=>u.*tan.(ϕ)/a)
       HoACos = diagm(0=>h./axcos)
       HxACos = diagm(0=>h.*axcos)
         dUdy = diagm(0=>Dy*u);

    # Form 1L-RSW Matrix
    #   [u1, v1, h1]

    A = [  UoACos                      (-F + dUdy - UxTanoA)[:, 2:N]     g*IoACos;
          -(F + 2*UxTanoA)[2:N,:]/kλ2                 UoACos[2:N,2:N]   -g/kλ2*Dy[2:N,:];
           HoACos                         (IoACos*Dy*HxACos)[:,  2:N]      UoACos] ;
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