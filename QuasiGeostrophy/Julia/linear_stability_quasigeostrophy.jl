using LinearAlgebra
using Plots
using IJulia

#Where to put these files?
include("cheb.jl")
include("parameters.jl")
include("geometry.jl")
include("profile.jl")
include("build_matrices.jl")
include("compute_spectrum.jl")
include("mesh.jl")
include("plot_fields.jl")

params=parameters(
     H = 1, 
    Ly = 20, 
    f₀ = 1, 
     g = 10, 
    Lj = 1, 
    Uj = 0.1,
    Ny = 256,
    dk = 5e-2,
    kₘ = 2,
    );

params.Fr = (params.f₀ * params.Lj)^2/(params.g * params.H);

# Wavenumbers
ks = collect(params.dk:params.dk:params.kₘ) / params.Lj;
Nk = length(ks);

# Add FD option
# Grid and basic state
         Ny = params.Ny;
Dy,  Dy2, y = geometry(Ny, params.Ly, "cheb");
dQdy,  U, Ψ = profile(  y, params);

plot_basic_state(y, dQdy, U, Ψ, params.Ly, "basic_state.png")

# initialize fields to store
Nmodes = 2;
     σ = zeros(Nmodes, Nk);
     ω = zeros(Nmodes, Nk);   
σmodes = zeros(ComplexF64, Ny+1, Nmodes, Nk);

# Do we want a structure with soln.U etc?
# Compute growth rates
for cnt in 1:length(ks)
    local k = ks[cnt]
    local A, B = build_AB(k, dQdy, U, Ψ, Dy, Dy2, y, params);
    local σ[:,cnt], ω[:,cnt], σmodes[2:end-1,:,cnt] = compute_spectrum(A, B, k, params, Nmodes);
end

plot_growth_rates(ks, σ, Nmodes, "growth_rates.png")

mode_number = 1;
      σ_max = maximum(σ[mode_number,:]);
    k_index = sortperm(σ[mode_number,:],rev=true)[1];
          k = ks[k_index]; # pick wavenumber

plot_1D_fields(k_index, k, Ny, y, σmodes, mode_number, "modes_1D.png")
plot_2D_fields(k_index, k, Ny, y, σmodes, mode_number, "modes_2d.png")