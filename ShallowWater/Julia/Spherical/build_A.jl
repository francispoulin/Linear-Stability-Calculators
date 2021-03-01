### Notation
#   small letters => vectors
# capital letters => matrices

function build_A(k, u, η, D, D2, ϕ, params)
    
    a = params.a
    N = params.Nϕ
    g = params.g
    
    k2 = k^2
    ik = 1im*k    

    h = η .+ params.H
    axcos = a * cos.(ϕ)
    
           D = D/a
       IxCos = diagm(0=>axcos)
       IoCos = diagm(0=>1 ./ axcos)
    TwoΩxSin = diagm(0=>params.TwoΩ*sin.(ϕ));
       UoCos = diagm(0=>u./axcos)
       UxTan = diagm(0=>u./axcos.*sin.(ϕ))
       HoCos = diagm(0=>h./axcos)
       HxCos = diagm(0=>h.*axcos)
          dU = diagm(0=>D*u);
          dH = diagm(0=>D*h);

    # Form 1L-RSW Matrix
    #   [u1, v1, h1]

     U = diagm(0 => u)
     H = diagm(0 => h)
     F = diagm(0 => params.TwoΩ .+ 0*ϕ)
    dU = diagm(0 => D*u)
    
    A = [                UoCos             (-TwoΩxSin - UxTan + dU)[:, 2:N]    g*IoCos;
          -(TwoΩxSin + 2*UxTan)[2:N,:]/k2                     UoCos[2:N,2:N]  -g/k2*D[2:N,:];
                         HoCos                  (IxCos * D * HxCos)[:, 2:N]     UoCos];
    return A
end

