mutable struct parameters
      H::Float64
     Ly::Float64
     f₀::Float64
      g::Float64
     Lj::Float64
     Uj::Float64
     Ny::Int64
     dk::Float64
     kₘ::Float64
end

function parameters(;
    H  = 500,              
    Ly = 2e5,             
    f₀ = 1e-4,            
    g  = 9.81,
    Lj = 1e4,              
    Uj = 0.1,
    Ny = 100,
    dk = 5e-6,
    kₘ = 2e-4
    )             
    return parameters(H, Ly, f₀ , g, Lj, Uj, Ny, dk, kₘ)
end