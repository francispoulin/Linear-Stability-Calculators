mutable struct parameters
      H::Float64            # meters               Physical
     Ly::Float64            # meters
     Lϕ::Float64            # radians
   ϕmin::Float64            # radials         
   ϕmid::Float64            # radians
   ϕmax::Float64            # radials         
      a::Float64            # meters

   TwoΩ::Float64            # 1/second             Coriolis
     f₀::Float64            # 1/second

      g::Float64            # meters/second^2      Buoyancy

      N::Int64              # Integer              Numerical

     Lj::Float64            # meters               Jet (specify profile)
     Uj::Float64            # meters/second

     dk::Float64            # 1/meters             Wavenumbers
   kmin::Float64            # 1/metres
   kmax::Float64            # 1/meters

     Fr::Float64            # Nondimensional       Nondimensional
     Ro::Float64            # Nondimensional
end

"""

Cartesian version of parameters function.

"""
function parameters(geometry::Cartesian;
     H = 1,
    Ly = 20,
    Lϕ = nothing,
  ϕmin = nothing,
  ϕmid = nothing,
  ϕmax = nothing,
     a = nothing,

  TwoΩ = nothing,
    f₀ = 1,

     g = 9.81,

     N = 100,

    Lj = 1,
    Uj = 0.1,

    dk = 1e-1,
  kmin = dk,
  kmax = 2,

    Fr = 0.0,
    Ro = 0.0
    )
    return parameters(H, Ly, f₀, g, N, Lj, Uj, N, dk, kmin, kmax, Fr, Ro)
end

"""

Spherical version of parameters function.

"""
function parameters(geometry::Spherical;
     H = 1,
    Ly = nothing,
    Lϕ =  π/6,
  ϕmin =  π/8,
  ϕmid =  π/4,
  ϕmax = 3π/8,
     a = 1,

  TwoΩ = 1,
    f₀ = nothing,

     g = 9.81,

     N = 100,

    Lj = 1,
    Uj = 0.1,

    dk = 1e-1,
  kmin = dk,
  kmax = 2,

    Fr = 0.0,
    Ro = 0.0
  )
  return parameters(H, Lϕ, ϕmin, ϕmid, ϕmax, a, TwoΩ, g, N, Lj, Uj, dk, kmin, kmax, Fr, Ro)
end
