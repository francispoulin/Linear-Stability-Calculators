function geometry(params)

    Ny = params.Ny
    Ly = params.Ly

    Dy,y = cheb(Ny); 
       y = y*(0.5*Ly);
      Dy = (2/Ly)*Dy;
     Dy2 = Dy*Dy;
    
    return Dy, Dy2, y
end
