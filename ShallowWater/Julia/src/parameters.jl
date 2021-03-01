mutable struct Parameters{I, F}
      H::F            # meters               Physical
     Ly::F            # meters
     Lϕ::F            # radians
   ϕmin::F            # radials         
   ϕmid::F            # radians
   ϕmax::F            # radials         
      a::F            # meters

   TwoΩ::F            # 1/second             Coriolis
     f₀::F            # 1/second

      g::F            # meters/second^2      Buoyancy

      N::I            # Integer              Numerical

     Lj::F            # meters               Jet (specify profile)
     Uj::F            # meters/second

     dk::F            # 1/meters             Wavenumbers
   kmin::F            # 1/metres
   kmax::F            # 1/meters

     Fr::F            # Nondimensional       Nondimensional
     Ro::F            # Nondimensional
end

"""

Cartesian version of parameters function.

"""
function Parameters(geometry::Cartesian;
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
    return Parameters(H, Ly, f₀, g, N, Lj, Uj, N, dk, kmin, kmax, Fr, Ro)
end

"""

Spherical version of parameters function.

"""
function Parameters(geometry::Spherical;
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
  return Parameters(H, Lϕ, ϕmin, ϕmid, ϕmax, a, TwoΩ, g, N, Lj, Uj, dk, kmin, kmax, Fr, Ro)
end
