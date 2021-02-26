using Cubature
using Statistics

function profile(ϕ, params);

    ϕmid = (ϕ[1] + ϕ[end])/2 
    Δη = params.Uj * params.Lj * params.TwoΩ / params.g

    U_function(x)  =  params.Uj * sech.((x .- ϕmid) / params.Lj).^2
    LHS(x) = - (params.TwoΩ .* sin.(x) .* params.a .+ U_function(x).*tan.(x)) .* U_function(x) / params.g
    
    U =  U_function(ϕ)
    E = zeros((params.Nϕ+1))
    for i in 1:params.Nϕ
        result = pquadrature(LHS, ϕ[end], ϕ[i])
        E[i] =  result[1]
    end
    E = E .- mean(E)

    return U,E
end




