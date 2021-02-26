function build_A(kx,U,E,Dy,Dy2,y,params)
    U1 = U
    h1 = E.+params.H

    # Define Differentiation Matrices for y: FD second order

    Ny  = length(U1)-1
    g   = params.g

    #  Define Matrices to set up EVP
    dU1 = diagm(0=>Dy*U1);
    dh1 = diagm(0=>Dy*h1);

    F    = diagm(0=>params.f₀.+0*y);
    #Id    = Matrix{Float64}(I, Ny+1, Ny+1);
    I    = diagm(0=>ones(Ny+1,1)[:])
    h1f  = diagm(0=>h1);

    kx2  = kx^2;
    ikx  = 1im*kx;
    
    # Form 1L-RSW Matrix
    #   [u1, v1, h1]

    A = [  diagm(0=>U1)  -F[:,2:Ny].+dU1[:,2:Ny]                  g*I;
         -F[2:Ny,:]/kx2       diagm(0=>U1[2:Ny])    -g/kx2*Dy[2:Ny,:];
           diagm(0=>h1)           Dy*h1f[:,2:Ny]         diagm(0=>U1)];
    
    return A
end