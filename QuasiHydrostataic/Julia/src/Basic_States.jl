using Cubature
using Statistics

Bickley(x; U, L, x₀) =  U * sech.((x .- x₀) / L).^2

mutable struct Bickley_Jet{V, P}
 phys::P
    u::V                            # meters/second
    η::V                            # meters
end

"""
Cartesian version.
"""
function Bickley_Jet(geometry::Cartesian;
    phys,
    Uj,
    Ljscale
    )

    Lj = phys.grid.Ly/Ljscale

    LHS(x) = - phys.f₀ / phys.g * Bickley(x; U=Uj, L=Lj, x₀=0)

    u = Bickley(phys.grid.y; U=Uj, L=Lj, x₀=0)
    η = zeros(size(phys.grid.y))
    for i in 1:phys.grid.N
        result = pquadrature(LHS, phys.grid.y[end], phys.grid.y[i])
        η[i] = result[1]
    end
    η = η .- mean(η)

    return Bickley_Jet(phys, u, η)
end

"""
Spherical version.
"""
function Bickley_Jet(geometry::Spherical;
    phys,
    Uj,
    Ljscale
    )

    Lj = phys.grid.Lϕ/Ljscale

    LHS(x) = - (phys.TwoΩ .* sin.(x) .* phys.grid.a .+ Bickley(x; U=Uj, L=Lj, x₀=phys.grid.ϕmid).*tan.(x)) .* Bickley(x; U=Uj, L=Lj, x₀=phys.grid.ϕmid) / phys.g

    u =  Bickley(phys.grid.ϕ; U=Uj, L=Lj, x₀ = phys.grid.ϕmid)
    η = zeros(size(phys.grid.ϕ))
    for i in 1:phys.grid.N
        result = pquadrature(LHS, phys.grid.ϕ[end], phys.grid.ϕ[i])
        η[i] =  result[1]
    end
    η = η .- mean(η)
   
    return Bickley_Jet(phys, u, η)
end

# Define basic state as we do in the Oceananigans library? 
# Pick a bunch of different profiles: U_shear, D_shear, ∫_shear