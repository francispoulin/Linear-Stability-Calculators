mutable struct parameters
      H::Float64
     Lϕ::Float64
   ϕmin::Float64
   TwoΩ::Float64
      g::Float64
      a::Float64
     Lj::Float64
     Uj::Float64
     Nϕ::Int64
     dk::Float64
     kₘ::Float64
     Fr::Float64
     Ro::Float64
end

function parameters(;
    H  = 500,              
    Lϕ = π/3,
  ϕmin = π/6,
  TwoΩ = 1e-4,            
    g  = 9.81,
    a  = 1,
    Lj = 1e4,              
    Uj = 0.1,
    Nϕ = 100,
    dk = 5e-6,
    kₘ = 2e-4,
    Fr = 0.0,   
    Ro = 0.0
    )             
    return parameters(H, Lϕ, ϕmin, TwoΩ, g, a, Lj, Uj, Nϕ, dk, kₘ, Fr, Ro)
end
