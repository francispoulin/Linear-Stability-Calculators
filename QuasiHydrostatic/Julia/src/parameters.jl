using Parameters

@with_kw struct Files{S}
  nc::S          = "data_spectrum.nc"
  json::S        = "data_parameters.json"
  plotgrowth::S   = "growth.png"
  plotslicem::S   = "growth_slice_m.png"
  plotmodes1D::S  = "plot_modes1D"
  plotmodes2D::S  = "plot_modes2D"
end

@with_kw struct Physics{F, I}
  Ω::F   = 2 * π / (24*3600)
  θ₀::F  = π/4
  N::F   = 1e-3
  g::F   = 9.81
  ν::F   = 0.26
  NT::I  = 1

  fz::F  = 2 * Ω * sin( θ₀ )
  fy::F  = 2 * Ω * cos( θ₀ ) * NT
end

include("cheb.jl")

mutable struct BuildMesh{F, I, V, M}
  L::  F
  N::  I 
  y::  V
  D::  M
  D2:: M
end

"""
Insert
"""
function BuildMesh(; 
    L = 1000e3, 
    N = 100
  )

  D, y = cheb(N)
  y    = L / 2 * y
  D    = 2 / L * D
  D2   = D * D
  
  return BuildMesh(L, N, y, D, D2)
end

mutable struct Grid{F, I, AM, M, V}
  Ly::F
  Lz::F
  Ny::I
  θ₀::F
  method::AM
  D:: M
  y:: V
end

"""
Insert
"""
function Grid(;
  Ly     = 1000e3,
  Lz     = 3e3,
  Ny     = 100,
  θ₀     = π/4,
  method = Cheb()
)

  M = BuildMesh(L = Ly, N = Ny)

return Grid(Ly, Lz, Ny, θ₀, method, M.D, M.y)
end

#function Bickley(grid.y, L, Umax, )

mutable struct Jet{G, P, F, S, V}
  grid::G         
  physics::P         
  L::F            # meters               Physical
  L_center::F     # meters
  H::F            # meters
  H_center::F     # meters
  Umax::F         # radians
  profile::S      
  
  ∫U::V
   U::V
  dU::V
end

"""
Insert
"""
function Jet(
  grid, 
  physics,
  L        = grid.Ly/10,
  L_center = 0.,
  H        = Inf,
  H_center = grid.Lz/2,
  Umax     = 14.6,
  profile  = "Bickley"
)
  ∫U =     L * Umax *   tanh.( grid.y / L )
   U =         Umax ./  cosh.( grid.y / L ).^2
  dU = 2 / L * Umax .*  tanh.( grid.y / L ) ./ cosh.( grid.y / L ).^2

return Jet(grid, physics, L, L_center, H, H_center, Umax, profile, ∫U, U, dU)
end

using JSON, JSON3

function Output_Parameters(grid, physics, jet, file)

  fz = physics.fz
  N  = physics.N
  g  = physics.g
  ν  = physics.ν

  U = jet.Umax
  L = jet.L

  # Non-dimensional parameters
  phi0  = U * fz * L
  Ro    = U / (fz * L)
  Bu    = (N / fz)^2
  delt  = phi0/(g * L)
  gamma = g^2 * ν / (fz*phi0^2)

  # FJP: Is this the output you really want???
  # Output parameters in terminal
  print("\n")
  print("Nondimensional Parameters\n")
  print("=========================\n")
  print("Ro      = ", Ro, "\n") 
  print("Bu      = ", Bu, "\n")
  print("phi0    = ", phi0, "\n")
  print("delta   = ", delt, "\n")
  print("gamma   = ", gamma, "\n")
  print("method  = ", grid.method, "\n")
  print("profile = ", jet.profile, "\n")
  print("\n")

  grid_dict    = Dict("Ly"=> grid.Ly, "Lz"=> grid.Lz, "Ny"=> grid.Ny, "θ₀"=>grid.θ₀)
  physics_dict = Dict("Ω"=> physics.Ω, "N"=> physics.N, "g"=> physics.g, "ν"=> physics.ν, "fz"=> physics.fz, "fy"=>physics.fy)
  jet_dict     = Dict("L"=> jet.L, "L_center"=> jet.L_center, "H"=> jet.H, "H_center"=> jet.H_center, "Umax"=> jet.Umax, "profile"=> jet.profile)
  big_dict     = Dict("grid"=> grid_dict, "physics"=> physics_dict, "jet"=>jet_dict)
  json_string  = JSON.json(big_dict)
  json_read    = JSON3.read(json_string)
  
  #JSON3.pretty(JSON3.write(json_read))
  
  open(file, "w") do io
    JSON3.pretty(io, json_read)
  end
  
end

using NCDatasets
using DataStructures

function save_spectrum(ωs, modes, ks, ms, Neigs, y, Ny, file)

  print("\n--> Saving the spectrum and modal structures in ", file, "\n")

  ys = [y; y[2:end-1]; y]

  ds = Dataset(file, "c")
  
  defDim(ds, "ks",    length(ks))
  defDim(ds, "ms",    length(ms))
  defDim(ds, "iEigs", Neigs)
  defDim(ds, "ys",    3*Ny+1)
  
  ds.attrib["title"] = "spectrum from linear stability analysis of inertial instability"

  # Coordinates
  k    = defVar(ds, "ks",    Float64, ("ks",),    attrib = OrderedDict("units" => "1/meter"))
  m    = defVar(ds, "ms",    Float64, ("ms",),    attrib = OrderedDict("units" => "1/meter"))
  iEig = defVar(ds, "iEigs", Float64, ("iEigs",), attrib = OrderedDict("units" => "None"))
  y    = defVar(ds, "ys",    Float64, ("ys",),    attrib = OrderedDict("units" => "m"))

  # Frequencies
  ωs_real = defVar(ds, "ω_real", Float64, ("iEigs", "ms", "ks"), attrib = OrderedDict("units" => "1/second"))
  ωs_real.attrib["comments"] = "frequencies that depend on k, m and iEigs"
  ωs_imag = defVar(ds, "ω_imag", Float64, ("iEigs", "ms", "ks"), attrib = OrderedDict("units" => "1/second"))
  ωs_imag.attrib["comments"] = "growth rates that depend on k, m and iEigs"

  # Modal Structures: ubv
  modes_real = defVar(ds, "uvb_real", Float64, ("ys", "iEigs", "ms", "ks"), attrib = OrderedDict("units" => "m/s"))
  modes_real.attrib["comments"] = "real part of modal structure that depend on k, m, ys and iEigs"
  modes_imag = defVar(ds, "uvb_imag", Float64, ("ys", "iEigs", "ms", "ks"), attrib = OrderedDict("units" => "m/s"))
  modes_imag.attrib["comments"] = "imaginary part of modal structure that depend on k, m, ys and iEigs"

  k[:]    = ks
  m[:]    = ms
  iEig[:] = collect(1:1:Neigs)
  y[:]    = ys

  ωs_real[:, :, :]       = real(ωs)
  ωs_imag[:, :, :]       = imag(ωs)
  modes_real[:, :, :, :] = real(modes)
  modes_imag[:, :, :, :] = imag(modes)

  close(ds)
end

function load_spectrum(file)

  ds = Dataset(file, "r")

  ωs    = ds["ω_real"][:, :, :]      + 1im*ds["ω_imag"][:, :, :]
  modes = ds["uvb_real"][:, :, :, :] + 1im*ds["uvb_imag"][:, :, :, :]

  Ny    = Int((size(modes)[1] - 1)/3)

  k     = ds["ks"][:]
  m     = ds["ms"][:]
  iEig  = Int.((ds["iEigs"][:]))
  y     = ds["ys"][1:Ny+1]
  
  close(ds)

  return ωs, modes, k, m, iEig, y
end

function read_parameters(file)

  json_string = read(file, String)

  JSON3.pretty(json_string)



end