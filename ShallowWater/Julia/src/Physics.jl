mutable struct Physics{F, V, G}
    grid::G
    TwoΩ::F            # 1/second             Coriolis
      f₀::F            # 1/second
       β::F            # 1/(second * meter)
       g::F            # meters/second^2      Buoyancy

coriolis::V            # vector with units 1/second

      Fr::F            # Nondimensional       Nondimensional
      Ro::F            # Nondimensional
end

"""

Cartesian version of Physics

"""
function Physics(geometry::Cartesian;
    grid,
    TwoΩ = 4π/(24*3600),
      f₀ = TwoΩ*sin(grid.ϕmid),
       β = TwoΩ*cos(grid.ϕmid)/grid.a,
       g = 9.81,
      Fr = 0.0,
      Ro = 0.0
)

    coriolis = f₀ .+ β * grid.y 

return Physics(grid, TwoΩ, f₀, β, g, coriolis, Fr, Ro)
end

"""

Spherical version of Physics

"""
function Physics(geometry::Spherical;
    grid, 
    TwoΩ = 4π/(24*3600),
      f₀ = TwoΩ*sin(grid.ϕmid),
       β = TwoΩ*cos(grid.ϕmid)/grid.a,
       g = 9.81,
      Fr = 0.0,
      Ro = 0.0
)

    coriolis = TwoΩ * sin.(grid.ϕ) 

return Physics(grid, TwoΩ, f₀, β, g, coriolis, Fr, Ro)
end

