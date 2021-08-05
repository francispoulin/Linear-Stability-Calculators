### Notation
#   small letters => vectors
# capital letters => matrices

function build_A(k, u, η, D, y, params)

    N = params.Ny
    g = params.g

    k2 = k^2
    ik = 1im*k

    # This removed the free-surface deformation in the basic state,
    # which is more like what is done in the Hydrostatic case
    
    h = 0*η .+ params.H

     U = diagm(0 => u)
     H = diagm(0 => h)
    dU = diagm(0 => D*u)
    dH = diagm(0 => D*h)
     F = diagm(0 => params.f₀ .+ 0*y)
     I = diagm(0 => ones(N+1))

    # Form 1L-RSW Matrix
    #   [u1, v1, h1]

    A = [    k*U         1im*(F-dU)[:, 2:N]    g*k*I;
          -1im*F[2:N,:]     k*U[2:N,2:N]    -1im*g*D[2:N,:];
             k*H         -1im*D*H[:,  2:N]       k*U];
    
    return A
end
