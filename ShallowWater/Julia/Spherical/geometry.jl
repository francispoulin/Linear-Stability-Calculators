function geometry(params)

       N = params.Nϕ
       L = params.Lϕ
    ϕmin = params.ϕmin
    
    D,y = cheb(N); 
      ϕ = L / 2 * (y .+ 1);
      D = 2 / L * D;
     D2 = D*D;
      ϕ = ϕ .+ ϕmin
    
    return D, D2, ϕ
end
