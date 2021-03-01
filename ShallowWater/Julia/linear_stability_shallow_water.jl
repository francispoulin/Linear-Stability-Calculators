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
include("Plotting_Functions.jl")

gridC = Grid(geometryC; N = 100)
physC = Physics(geometryC; grid=gridC, β=0.0)

gridS = Grid(geometryS; N = 100)
physS = Physics(geometryS; grid=gridS, β=0.0)

# Wavenumbers
  dk = 5e-6                         # 1/meters
kmax = 2e-4
  ks = collect(dk:dk:kmax)
  Nk = length(ks)

        Uj = 1.0                   # meters/second
   Ljscale = 20                  # meters
backgroundC = Bickley_Jet(geometryC; phys=physC, Uj, Ljscale)
backgroundS = Bickley_Jet(geometryS; phys=physS, Uj, Ljscale)

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

include("src/Compute_Eigenvalues.jl")

for (index, k) in enumerate(ks)
     print("index = $index and k = $k\n")
     AC = Build_Matrix(k=k, solution=backgroundC, phys=physC)
end

#=
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

#string.(fieldnames(typeof(gridC)))
