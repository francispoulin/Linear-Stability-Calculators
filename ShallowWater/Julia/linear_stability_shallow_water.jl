"""
Linear-Stability-Calculators: Shallow Water Model
=================================================
"""

include("src/shallow_water_setup.jl")

             ks = Dict()
ks[Cartesian()] = collect(2e-5:2e-5:2.5e-3)
ks[Spherical()] = collect(2e-5:2e-5:2.5e-3) * 6.3781e6 * cos(π/4)

        Uj = 1.0                 # Jet Parameters
   Ljscale = 20                 
    Nmodes = 2                    
         σ = Dict()
         ω = Dict()
    σmodes = Dict()
      phys = Dict()
background = Dict()

geometrys = (Cartesian(), Spherical()) 

for geometry in geometrys

     @info string("Setting Geometry = ", geometry)

                    grid = Grid(geometry, N=250)
          phys[geometry] = Physics(geometry; grid=grid)
     background[geometry] = Bickley_Jet(geometry; phys=phys[geometry], Uj, Ljscale)
  
     plot_basic_state(
          grid, 
          background, 
          geometry, 
          string("basic_state_",typeof(geometry),".png")
          )

     for (index, k) in enumerate(ks[geometry])
          A = build_matrix(geometry; k, background, phys)
          σ[(k,geometry)], ω[(k,geometry)], σmodes[(k,geometry)] = compute_spectrum(A, k, phys, Nmodes);
          @info string("     k = ", k, " and max growth = ", maximum(σ[(k,geometry)]), "\n")
     end

     plot_growth_rates(
          ks, 
          σ, 
          phys,
          geometry, 
          Nmodes, 
          string("growth_rates_",typeof(geometry),"_Bickley_Jet.png")
          )
     
          mode_number = 1;
          σ_max = maximum( σ[(k,geometry)][mode_number] for k in ks[geometry]);
        k_index = sortperm([σ[(k,geometry)][mode_number] for k in ks[geometry]],rev=true)[1];
              k = ks[geometry][k_index];

          print("σ_max = ",σ_max, " with k = ", k, " at k_index =", k_index, "\n")
          plot_1D_streamfunction(
               k_index, 
               k, 
               phys, 
               σmodes, 
               geometry, 
               mode_number, 
               string("modes_1D_streamfunction_",typeof(geometry),".png")
               )

          plot_2D_streamfunction(
               k_index, 
               k, 
               phys, 
               σmodes, 
               geometry,
               mode_number, 
               string("modes_2D_streamfunction_",typeof(geometry),".png")
               )

          plot_2D_vorticity(
               phys, 
               k_index, 
               k, 
               σmodes, 
               geometry,
               mode_number, 
               string("modes_2D_vorticity_",typeof(geometry),".png")
               )

end

#=
plot_1D_streamfunction(k_index, k, y, σmodes, mode_number, "modes_1D_streamfunction.png")
plot_2D_streamfunction(k_index, k, y, σmodes, mode_number, "modes_2D_streamfunction.png")
plot_2D_vorticity( Dy, k_index, k, y, σmodes, mode_number, "modes_2D_vorticity.png")
=#

#string.(fieldnames(typeof(gridC)))
