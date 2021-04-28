mutable struct Physics{F, V, G}
     grid::G
     TwoΩ::F            # 1/second             Coriolis
       fz::F            # 1/second
       fy::F            # 1/second   
       βz::F            # 1/(second * meter)
       βy::F            # 1/(second * meter)
        g::F            # meters/second^2      Buoyancy

coriolisz::V            # vector with units 1/second
coriolisy::V            # vector with units 1/second

       Bu::F            # Burger number 
       Ro::F            # Rossby number 
        γ::F            # Non-traditional parameter 
       Ek::F            # vertical Ekmann number 
end

"""

Cartesian version of Physics

"""
function Physics(geometry::Cartesian;
    grid,
    TwoΩ = 4π/(24*3600),
      fz = TwoΩ*sin(grid.ϕmid),
      fy = TwoΩ*cos(grid.ϕmid),
      βz = TwoΩ*cos(grid.ϕmid)/grid.a,
      βy = 2*TwoΩ*sin(grid.ϕmid)/grid.a,
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

# Add background stratification N²?  
# Bu = (N/f)²