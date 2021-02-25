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

# Wavenumbers
ks = collect(params.dk:params.dk:params.kₘ) / params.Lj;
Nk = length(ks);

# Grid and basic state
Dy, Dy2, y = geometry(params.Ny, params.Ly, "cheb");
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
    local A = build_A(k, U, E, Dy, Dy2, y, params);
    local σ[:,cnt], ω[:,cnt], σmodes[:,:,cnt] = compute_spectrum(A, k, params, Nmodes);
end

plot_growth_rates(ks, σ, Nmodes, "growth_rates_Cartesian_Bickley_jet.png")

mode_number = 1;
      σ_max = maximum(σ[mode_number,:]);
    k_index = sortperm(σ[mode_number,:],rev=true)[1];
          k = ks[k_index]; # pick wavenumber

plot_1D_streamfunction(k_index, k, params.Ny, y, σmodes, mode_number, "modes_1D_streamfunction.png")
plot_2D_streamfunction(k_index, k, params.Ny, y, σmodes, mode_number, "modes_2D_streamfunction.png")
plot_2D_vorticity(k_index, k, params.Ny, y, Dy2, σmodes, mode_number, "modes_2D_vorticity.png")

