using Cubature
using Statistics

# Hard coded for a Bickley jet but want to generalize to other profiles including
# -> Bump and Guassian jets
# -> Vortices of different profiles
# -> In Cartesian and Spherical profiles

# Velocity is specified to be a Bickley jet
# Free-Surface is computed using numerical integration

function profile(y, params);

    Δη = params.Uj * params.Lj * params.f₀ / params.g

    U_function(x)  =  params.Uj * sech.(x / params.Lj).^2

    U =  U_function(y)
    E = zeros((params.Ny+1))
    for i in 1:params.Ny
        result = pquadrature(U_function, y[end], y[i])
        E[i] = - params.f₀ / params.g * result[1]
    end
    E = E .- mean(E)

    return U,E
end
