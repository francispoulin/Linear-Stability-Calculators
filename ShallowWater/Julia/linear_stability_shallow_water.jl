"""

Linear-Stability-Calculators
============================
Repo:      Linear-Stability-Calculators
Code:      ShallowWater/Julia/linear_stability_shallow_water.jl
Model:     Rotating Shallow Water
Geometry:  Cartesian/Spherical
Structure: Bickley Jet

Nondimensional Parameters
=========================
Fr = Froude Number
Ro = Rossby Number

Created Feb 26, 2021
By Francis J. Poulin

"""

include("src/shallow_water_setup.jl")

geometryC = Cartesian()
geometryS = Spherical()

include("src/Grid.jl")
include("src/Physics.jl")
include("src/Basic_States.jl")

gridC = Grid(geometryC; N = 100)
physC = Physics(geometryC; gridC, β=0.0)

gridS = Grid(geometryS; N = 100)
physS = Physics(geometryS; gridS, β=0.0)

#=

# Wavenumbers
  dk = 5e-6                         # 1/meters
kmax = 2e-4
  ks = collect(dk:dk:kmax)
  Nk = length(ks)

        Uj = 1.0                   # meters/second
   Ljscale = 20                  # meters
backgroundC = Bickley_Jet(geometryC; gridC, physC, Uj, Ljscale)
backgroundS = Bickley_Jet(geometryS; gridS, physS, Uj, Ljscale)

fileC = string("basic_state_",typeof(geometryC),".png")
plot_basic_state(gridC.y, backgroundC, fileC)

fileS = string("basic_state_",typeof(geometryS),".png")
plot_basic_state(gridS.ϕ, backgroundS, fileS)

# initialize fields to store
Nmodes = 2
     σC = zeros(Nmodes, Nk);
     ωC = zeros(Nmodes, Nk);   
σmodesC = zeros(ComplexF64, 3*gridC.N+1, Nmodes, Nk);

     σS = zeros(Nmodes, Nk);
     ωS = zeros(Nmodes, Nk);   
σmodesS = zeros(ComplexF64, 3*gridS.N+1, Nmodes, Nk);
=#

#include("src/parameters.jl")

#=
include("parameters.jl")
include("geometry.jl")
include("profile.jl")
include("build_A.jl")
include("compute_spectrum.jl")
include("mesh.jl")
include("plot_fields.jl")
 
params=parameters(
     H = 1, 
    Ly = 20, 
    f₀ = 1, 
     g = 10, 
    Lj = 1, 
    Uj = 1.0,
    Ny = 250,
    dk = 5e-2,
    kₘ = 2
    );
    
params.Fr = (params.f₀ * params.Lj)^2/(params.g * params.H);
params.Ro =  params.Uj/(params.f₀ * params.Lj)

print("\n")
print("Linear-Stability-Caluculator\n")
print("============================\n\n")
print("Repo:      Linear-Stability-Calculators\n")
print("Code:      ShallowWater/Julia/linear_stability_shallow_water.jl\n")
print("Model:     Rotating Shallow Water\n")
print("Geometry:  Cartesian\n")
print("Structure: Bickley Jet\n")
print("\n")
print("Nondimensional Parameters\n")
print("=========================\n")
print("Fr = ", params.Fr, "\n")
print("Ro = ", params.Ro, "\n")

# Wavenumbers
ks = collect(params.dk:params.dk:params.kₘ) / params.Lj;
Nk = length(ks);

# Grid and basic state
Dy, Dy2, y = geometry(params);
U,  E      = profile(         y, params);

plot_basic_state(y, U, E, params.Ly, "basic_state.png")

# initialize fields to store 
Nmodes = 2
     σ = zeros(Nmodes, Nk);
     ω = zeros(Nmodes, Nk);   
σmodes = zeros(ComplexF64, 3*params.Ny+1, Nmodes, Nk);

# Compute growth rates
for cnt in 1:Nk
    local k = ks[cnt]
    local A = build_A(k, U, E, Dy, y, params);
    local σ[:,cnt], ω[:,cnt], σmodes[:,:,cnt] = compute_spectrum(A, k, params, Nmodes);
end

plot_growth_rates(ks, σ, Nmodes, "growth_rates_Cartesian_Bickley_jet.png")

mode_number = 1;
      σ_max = maximum(σ[mode_number,:]);
    k_index = sortperm(σ[mode_number,:],rev=true)[1];
          k = ks[k_index]; # pick wavenumber

plot_1D_streamfunction(k_index, k, y, σmodes, mode_number, "modes_1D_streamfunction.png")
plot_2D_streamfunction(k_index, k, y, σmodes, mode_number, "modes_2D_streamfunction.png")
plot_2D_vorticity( Dy, k_index, k, y, σmodes, mode_number, "modes_2D_vorticity.png")
=#