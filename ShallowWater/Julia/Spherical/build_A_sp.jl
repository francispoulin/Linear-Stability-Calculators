function build_A_sp(kx,U,E,Dy,Dy2,y,parms)
    U1 = U
    h1 = E.+parms.H

    # Define Differentiation Matrices for y: FD second order

    Ny  = length(U1)-1
    g1  = parms.g1

    #  Define Matrices to set up EVP
  
    dU1 = Diagonal(Dy*U1);
    dh1 = Diagonal(Dy*h1);

    F    = Diagonal(parms.f0.+0*y);
    Id    = Matrix{Float64}(I, Ny+1, Ny+1);
    h1f  = Diagonal(h1);

    kx2  = kx^2;
    ikx  = 1im*kx;


    # Form 1L-RSW Matrix
    #   [u1, v1, h1]

    A = [[Diagonal(U1),      -F[:,2:Ny].+dU1[:,2:Ny],                g1*I], 
        [-F[2:Ny,:]/kx2,         Diagonal(U1[2:Ny]),  -g1/kx2*Dy[2:Ny,:]],
        [Diagonal(h1),       Dy*h1f[:,2:Ny],                    Diagonal(U1)]];
    return sparse(A)
end