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

# Wavenumbers: Put in a dictionary?
dk = 2e-5                         # 1/meters
kmax = 2e-4
ks = collect(dk:dk:kmax)
Nk = length(ks)
#ksS = ksC * gridC.a /1000

Uj = 1.0                      # Jet Parameters
Ljscale = 20                 
Nmodes = 2                    # Number of modes to keep

geometrys = (Cartesian(), Spherical())
for geometry in geometrys

     grid = Grid(geometry; N = 100)
     phys = Physics(geometry; grid=grid, β=0.0)
     background = Bickley_Jet(geometry; phys=phys, Uj, Ljscale)     # turn into a function?
  
     file = string("basic_state_",typeof(geometry),".png")
     plot_basic_state(grid.y, background, file)

     σ = zeros(Nmodes, Nk);     # Dictionary!!
     ω = zeros(Nmodes, Nk);   
     σmodes = zeros(ComplexF64, 3*grid.N+1, Nmodes, Nk);

     # Put in a parallel for loop and Use dictionarys
     for (index, k) in enumerate(ks)
          print("index = $index and k = $k\n")
          A = build_matrix(geometry; k=k, solution=background, phys=phys)
          σ[:,index], ω[:,index], σmodes[:,:,index] = compute_spectrum(A, k, phys, Nmodes);
     end

     file = string("growth_rates_",typeof(geometry),"_Bickley_Jet.png")
     plot_growth_rates(ks, σ, Nmodes, file)
            
end

#=
mode_number = 1;
      σ_max = maximum(σ[mode_number,:]);
    k_index = sortperm(σ[mode_number,:],rev=true)[1];
          k = ks[k_index]; # pick wavenumber

plot_1D_streamfunction(k_index, k, y, σmodes, mode_number, "modes_1D_streamfunction.png")
plot_2D_streamfunction(k_index, k, y, σmodes, mode_number, "modes_2D_streamfunction.png")
plot_2D_vorticity( Dy, k_index, k, y, σmodes, mode_number, "modes_2D_vorticity.png")
=#

#string.(fieldnames(typeof(gridC)))
