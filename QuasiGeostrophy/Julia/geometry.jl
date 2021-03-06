function geometry(Ny, Ly, method)

    if method=="cheb"
        Dy,y   = cheb(Ny); 
        y      = y*(0.5*Ly);
        Dy     = (2/Ly)*Dy;
        Dy2    = Dy*Dy;
        return Dy, Dy2, y
    end

    #=
    if method=="FD2" # BROKEN
        y = range(-Ly/2, Ly/2, length=Ny+1);
        hy=y[2]-y[1];
        e = ones(Ny+1,1);
        Dy = Tridiagonal(-1*e[1:Ny]/(2*hy), 0*e[:], e[1:Ny]/(2*hy));
        Dy[1,1:2] = [-1 1]/hy;
        Dy[Ny+1,Ny:Ny+1] = [-1 1]/hy;
        Dy2 = Tridiagonal(e[1:Ny]/(hy^2), -2*e[:]/(hy^2), e[1:Ny]/(hy^2));
        #Dy2[1,1:3] = [1 -2 1]/hy^2;
        #Dy2[Ny+1,Ny-1:Ny+1] = [1 -2 1]/hy^2;  
        return Dy,Dy2,y 
    end 
    =#

end