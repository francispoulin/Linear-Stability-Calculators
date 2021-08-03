using LinearAlgebra
using Plots

include("cheb.jl")
include("parameters.jl")
include("geometry.jl")
include("profile.jl")
include("build_A.jl")
include("compute_spectrum.jl")
include("mesh.jl")
include("plot_fields.jl")

Omega = 2*π/(24*3600)
lat   = π/4

params=parameters(
     H = 2e3, 
    Ly = 600e3, 
    f₀ = 2*Omega*cos(lat), 
     g = 10, 
    Lj = 60e3, 
    Uj = 1.0,
    Ny = 200,
    dk = 5e-2,
    kₘ = 2e0
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
