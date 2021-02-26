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
 
params=parameters(
     H = 1, 
    Lϕ = π/100,
  ϕmin = π/4,
   TwoΩ = 1, 
      g = 10,
      a = 1,
    Lj = π/2000, 
    Uj = 1.0,
    Nϕ = 100,
    dk = 1e-2,
    kₘ = 3
    ); 
    
params.Fr = (params.TwoΩ * params.Lj)^2/(params.g * params.H);
params.Ro =  params.Uj/(params.TwoΩ * params.Lj)

print("\n")
print("Linear-Stability-Caluculator\n")
print("============================\n\n")
print("Repo:      Linear-Stability-Calculators\n")
print("Code:      ShallowWater/Julia/Spherical/linear_stability_shallow_water.jl\n")
print("Model:     Rotating Shallow Water\n")
print("Geometry:  Spherical\n")
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
Dϕ, Dϕ2, ϕ = geometry(params);
      U,  E = profile(  ϕ, params);

plot_basic_state(ϕ, U, E, "basic_state.png")

# initialize fields to store
Nmodes = 2
     σ = zeros(Nmodes, Nk);
     ω = zeros(Nmodes, Nk);   
σmodes = zeros(ComplexF64, 3*params.Nϕ+1, Nmodes, Nk);

# Compute growth rates
for cnt in 1:Nk
    local k = ks[cnt]
    local A = build_A(k, U, E, Dϕ, Dϕ2, ϕ, params);
    local σ[:,cnt], ω[:,cnt], σmodes[:,:,cnt] = compute_spectrum(A, k, params, Nmodes);
end

plot_growth_rates(ks, σ, Nmodes, "growth_rates_spherical_Bickley_jet.png")

#print("σ = ", σ[1,:],"\n")

mode_number = 1;
      σ_max = maximum(σ[mode_number,:]);
    k_index = sortperm(σ[mode_number,:],rev=true)[1];
          k = ks[k_index]; # pick wavenumber

plot_1D_streamfunction(k_index, k, ϕ, σmodes, mode_number, "modes_1D_streamfunction.png")
plot_2D_streamfunction(k_index, k, ϕ, σmodes, mode_number, "modes_2D_streamfunction.png")
plot_2D_vorticity(Dϕ,  k_index, k, ϕ, σmodes, mode_number, "modes_2D_vorticity.png")
