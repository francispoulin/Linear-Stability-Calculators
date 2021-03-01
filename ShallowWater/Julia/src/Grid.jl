include("cheb.jl")

mutable struct Grid{I, F, V, M}
    H::F            # meters               Physical
   Ly::F            # meters
   Lϕ::F            # radians
 ϕmid::F            # radians
    a::F            # meters
    N::I            # Integer              Numerical
    y::V 
    ϕ::V 
    D::M
   D2::M
end

"""

Cartesian version of Grid function.

"""
function Grid(geometry::Cartesian;
     H = 500.0,
    Ly = 200e3,
    Lϕ = 3.135731e-2,
  ϕmid = π/4,
     a = 6.3781e6,
     N = 100
 )

 D, y = cheb(N)
    y = Ly / 2 * y
    ϕ = ϕmid .+ y/a
    D = 2 / Ly * D
   D2 = D * D

 return Grid(H, Ly, Lϕ, ϕmid, a, N, y, ϕ, D, D2)
end

"""

Spherical version of Grid function.

"""
function Grid(geometry::Spherical;
     H = 500.0,
    Ly = 200e3,
    Lϕ = 3.135731e-2,
  ϕmid = π/4,
     a = 6.3781e6,
     N = 100
)

D, ϕ = cheb(N)
   ϕ = Lϕ / 2 * ϕ
   y = ϕ * a
   ϕ = ϕ .+ ϕmid
   D = 2 / Lϕ * D
  D2 = D * D

return Grid(H, Ly, Lϕ, ϕmid, a, N, y, ϕ, D, D2)
end

