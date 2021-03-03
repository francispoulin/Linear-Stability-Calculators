using LinearAlgebra
using Plots
using IJulia
using Printf

function meshgrid(x, y)
    x = reshape(x, 1, length(x))
    y = reshape(y, length(y), 1)
    X = repeat(x, length(y), 1)
    Y = repeat(y, 1, length(x))
    return X, Y
end

# choose y label to be either y or ϕ, depending on the geometry
function plot_basic_state(grid, background, geometry, file)

    if geometry == Cartesian()
        y = grid.y/1e3
        xlabel = "y (km)"
    elseif geometry == Spherical()
        y = rad2deg.(grid.ϕ)
        xlabel = "ϕ (°)"
    end 

    solution = background[geometry]

     p1 = plot( 
         y, 
         solution.u,
         linewidth = 3,
         label = "u",
         title = "Velocity",
         ylabel = "(m/s)",
         xlims = (y[end], y[1])
         )
         
     p2 = plot(
         y, 
         solution.η,
         linewidth = 3,
         label = "η",
         xlabel = xlabel,
         ylabel = "(m)",
         title = "Free-Surface",
         xlims = (y[end], y[1])
         )
 
     plt1 = plot(p1, p2, layout = (2,1))    
     savefig(plt1, file)
end

function plot_growth_rates(ks, σ, phys, geometry, Nmodes, file)

    if geometry == Cartesian()
        k = ks[geometry]*1e3
        xlabel = "k (1/km)"
    elseif geometry == Spherical()
        k = ks[geometry]/1e3
        xlabel = "k (10^3/rad)"
    end 

    titles = ["1st mode", "2nd mode"]
    plt = plot()
    for (index,value) in enumerate(1:Nmodes)
        σnormalized = [σ[(k, geometry)][index] for k in ks[geometry]]/phys[geometry].f₀
        plot!(plt,
        k, 
        σnormalized,
        linewidth = 3,
        label=titles[index],
        xlabel=xlabel,
        ylabel="σ/f₀",
        title="Growth Rates",
        xlims=(k[1], k[end])
        )
    end
    savefig(plt,file)

end


#ticks = collect(0:0.2:1)
#ticklabels = [ @sprintf("%5.1f",x) for x in ticks ]
#plot(x,y)
#plot!(xticks=(ticks,ticklabels))

function  make_tick_labels(fs, Nfs)
    fs_max = maximum([real(fs); imag(fs)])
    fs_min = minimum([real(fs); imag(fs)])
    dfs = (fs_max - fs_min)/Nfs
  
    ticks = collect(fs_min:dfs:fs_max)
    ticklabels = [ @sprintf("%f",f) for f in ticks]

    return ticks, ticklabels
end

function plot_1D_streamfunction(k_index, k, phys, σmodes, geometry, mode_number, file)

    if geometry == Cartesian()
        y = phys[geometry].grid.y/1e3
        xlabel = "y (km)"
    elseif geometry == Spherical()
        y = rad2deg.(phys[geometry].grid.ϕ)
        xlabel = "ϕ (°)"
    end 

    N = phys[geometry].grid.N
    Lx = 2*pi/k;
    Nx = N;
     x = LinRange(0, Lx, Nx+1)
  X, Y = meshgrid(x,y)
  
  uvec = σmodes[(k,geometry)][    1:N+1,   mode_number];
  vvec = σmodes[(k,geometry)][  N+2:2*N+2, mode_number];
  ηvec = σmodes[(k,geometry)][2*N+1:3*N+1, mode_number];
  
    p1 = plot(
      y, 
      real(uvec),
      linewidth = 3, 
      label = "real",
      title = "u",
      yticks=make_tick_labels(uvec, 3),
      xlims = (y[end], y[1])
      )
  plot!(
      p1,
      y,
      linewidth = 3,
      imag(uvec),
      label = "imag")
  p2 = plot( 
      y, 
      real(vvec),
      linewidth = 3, 
      label = "real",
      title = "v",
      yticks=make_tick_labels(vvec, 3),
       xlims = (y[end], y[1])
      )
  plot!(
      p2,
      y,
      linewidth = 3,
      imag(vvec),
      label = "imag")
  p3 = plot(
      y, 
      real(ηvec),
      linewidth = 3, 
      label = "real",
      title = "η",
      yticks=make_tick_labels(ηvec, 3),
      xlabel = xlabel,
      xlims = (y[end], y[1])
      )
  plot!(
      p3,
      y,
      linewidth = 3,
      imag(ηvec),
      label = "imag")
  plt = plot(p1, p2, p3, layout = (3,1))   
  savefig(plt, file)
  #display(plt)
  
end


function plot_2D_streamfunction(k_index, k, phys, σmodes, geometry, mode_number, file)

    N = phys[geometry].grid.N
    Nx = N;
    Lx = 2*pi/k;
    x = LinRange(0, Lx, Nx+1)

    if geometry == Cartesian()
        y = phys[geometry].grid.y/1e3
        ylabel = "y (km)"
        x = x./1e3
        kscaled = k*1e3
        xlabel = "x (km)"
    elseif geometry == Spherical()
        y = rad2deg.(phys[geometry].grid.ϕ)
        ylabel = "ϕ (°)"
        x = rad2deg.(x)
        kscaled = k * π / 180
        xlabel = "λ (°)"
     end 

  X, Y = meshgrid(x,y)
  
  uvec = σmodes[(k,geometry)][    1:N+1,   mode_number];
  vvec = σmodes[(k,geometry)][  N+2:2*N+2, mode_number];
  ηvec = σmodes[(k,geometry)][2*N+1:3*N+1, mode_number];

    u =     repeat(real(uvec),1,Nx+1).*cos.(kscaled*X) -    repeat(imag(uvec),1,Nx+1).*sin.(kscaled*X);
    v = -k.*repeat(imag(vvec),1,Nx+1).*cos.(kscaled*X) - k.*repeat(real(vvec),1,Nx+1).*sin.(kscaled*X);
    η =     repeat(real(ηvec),1,Nx+1).*cos.(kscaled*X) -    repeat(imag(ηvec),1,Nx+1).*sin.(kscaled*X);

    kwargs = (
        xlabel = "x (km)",
        ylabel = ylabel,
          fill = true,
        levels = 20,
     linewidth = 0,
         color = :balance,
      colorbar = true,
          xlim = (x[1],   x[end]),
          ylim = (y[end], y[1])
    )

    #print("x ticks = ", make_tick_labels(x, 2), "\n")

    u_plt = contour(x, y[end:-1:1], u[:, end:-1:1], title="u"; kwargs...)
    v_plt = contour(x, y[end:-1:1], v[:, end:-1:1], title="v"; kwargs...)
    η_plt = contour(x, y[end:-1:1], η[:, end:-1:1], title="η"; kwargs...)

    plt1 = plot(u_plt, v_plt, η_plt, layout = (1,3))
    savefig(plt1, file)

end


function plot_2D_vorticity(phys, k_index, k, σmodes, geometry, mode_number, file)

    N = phys[geometry].grid.N
    Lx = 2*pi/k;
    Nx = N;
    x = LinRange(0, Lx, Nx+1)

    if geometry == Cartesian()
        y = phys[geometry].grid.y/1e3
        ylabel = "y (km)"
        Dy = phys[geometry].grid.D
        x = x./1e3
        kscaled = k*1e3
        xlabel = "x (km)"
    elseif geometry == Spherical()
        y = rad2deg.(phys[geometry].grid.ϕ)
        ylabel = "ϕ (°)"
        Dy = phys[geometry].grid.D/phys[geometry].grid.a
        x = rad2deg.(x)
        kscaled = k * π / 180
        xlabel = "λ (°)"
    end 

    X, Y = meshgrid(x,y)
  
  uvec = σmodes[(k,geometry)][    1:N+1,   mode_number];
  vvec = σmodes[(k,geometry)][  N+2:2*N+2, mode_number];

    ζvec = im * k * vvec - Dy * uvec    
    ζ    = repeat(real(ζvec),1,Nx+1).*cos.(kscaled*X) - repeat(imag(ζvec),1,Nx+1).*sin.(kscaled*X);

    kwargs = (
        xlabel = xlabel,
        ylabel = ylabel,
          fill = true,
        levels = 20,
     linewidth = 0,
         color = :balance,
      colorbar = true,
          xlim = (  x[1], x[end]),
         ylim = (y[end],  y[1])
    )

    plt = contour(x, y[end:-1:1], ζ[:, end:-1:1], title="ζ"; kwargs...)

    savefig(plt, file)

end

#string.(fieldnames(typeof(gridC)))
