using Cubature
using Statistics

Bickley(x; U, L, x₀) =  U * sech.((x .- x₀) / L).^2

mutable struct Bickley_Jet{V}
    u::V                            # meters/second
    η::V                            # meters
end

"""
Cartesian version.
"""
function Bickley_Jet(geometry::Cartesian;
    grid,
    phys,
    Uj,
    Ljscale
    )

    Lj = grid.Ly/Ljscale

    LHS(x) = - phys.f₀ / phys.g * Bickley(x; U=Uj, L=Lj, x₀=0)

    u = Bickley(grid.y; U=Uj, L=Lj, x₀=0)
    η = zeros(size(grid.y))
    for i in 1:grid.N
        result = pquadrature(LHS, grid.y[end], grid.y[i])
        η[i] = result[1]
    end
    η = η .- mean(η)

    return Bickley_Jet(u, η)
end

"""
Spherical version.
"""
function Bickley_Jet(geometry::Spherical;
    grid,
    phys,
    Uj,
    Ljscale
    )

    Lj = grid.Lϕ/Ljscale

    LHS(x) = - (phys.TwoΩ .* sin.(x) .* grid.a .+ Bickley(x; U=Uj, L=Lj, x₀=grid.ϕmid).*tan.(x)) .* Bickley(x; U=Uj, L=Lj, x₀=grid.ϕmid) / phys.g

    u =  Bickley_Profile(grid.ϕ; U=Uj, L=Lj, x₀ = grid.ϕmid)
    η = zeros(size(grid.ϕ))
    for i in 1:grid.N
        result = pquadrature(LHS, grid.ϕ[end], grid.ϕ[i])
        η[i] =  result[1]
    end
    η = η .- mean(η)
   
    return Bickley_Jet(u, η)
end