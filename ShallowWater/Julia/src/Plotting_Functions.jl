using Plots
using IJulia

using LinearAlgebra

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
        y = grid.y
        xlabel = "y (m)"
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
         xlims = (y[end], y[1])
         )
         
     p2 = plot(
         y, 
         solution.η,
         linewidth = 3,
         label = "η",
         xlabel = xlabel,
         title = "Free-Surface",
         xlims = (y[end], y[1])
         )
 
     plt1 = plot(p1, p2, layout = (2,1))    
     savefig(plt1, file)
end

function plot_growth_rates(ks, σ, geometry, Nmodes, file)

    titles = ["1st mode", "2nd mode"]
    plt = plot()
    for (index,value) in enumerate(1:Nmodes)
        plot!(plt,
        ks[geometry], 
        [σ[(k, geometry)][index] for k in ks[geometry]], 
        linewidth = 3,
        label=titles[index],
        title="Growth Rates",
        xlabel="k",
        xlims=(ks[geometry][1], ks[geometry][end])
        )
    end
    savefig(plt,file)

end



function plot_1D_streamfunction(k_index, k, phys, σmodes, geometry, mode_number, file)

    if geometry == Cartesian()
        y = phys[geometry].grid.y
        xlabel = "y (m)"
    elseif geometry == Spherical()
        y = rad2deg.(grid.ϕ)
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
      xlabel = xlabel,
      title = "u",
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
      xlims = (y[end], y[1])
      )
  plot!(
      p3,
      y,
      linewidth = 3,
      imag(vvec),
      label = "imag")
  plt = plot(p1, p2, p3, layout = (3,1))   
  savefig(plt, file)
  #display(plt)
  
end


function plot_2D_streamfunction(k_index, k, phys, σmodes, geometry, mode_number, file)

    if geometry == Cartesian()
        y = phys[geometry].grid.y
        ylabel = "y (m)"
    elseif geometry == Spherical()
        y = rad2deg.(grid.ϕ)
        ylabel = "ϕ (°)"
    end 

    N = phys[geometry].grid.N
    Nx = N;
    Lx = 2*pi/k;
    x = LinRange(0, Lx, Nx+1)
  X, Y = meshgrid(x,y)
  
  uvec = σmodes[(k,geometry)][    1:N+1,   mode_number];
  vvec = σmodes[(k,geometry)][  N+2:2*N+2, mode_number];
  ηvec = σmodes[(k,geometry)][2*N+1:3*N+1, mode_number];

    u =     repeat(real(uvec),1,Nx+1).*cos.(k*X) -    repeat(imag(uvec),1,Nx+1).*sin.(k*X);
    v = -k.*repeat(imag(vvec),1,Nx+1).*cos.(k*X) - k.*repeat(real(vvec),1,Nx+1).*sin.(k*X);
    η =     repeat(real(ηvec),1,Nx+1).*cos.(k*X) -    repeat(imag(ηvec),1,Nx+1).*sin.(k*X);

    kwargs = (
        xlabel = "x",
        ylabel = ylabel,
          fill = true,
        levels = 20,
     linewidth = 0,
         color = :balance,
      colorbar = true,
          xlim = (    0, Lx),
         ylim = (y[end], y[1])
    )

    u_plt = contour(x, y[end:-1:1], u[:, end:-1:1], title="u"; kwargs...)
    v_plt = contour(x, y[end:-1:1], v[:, end:-1:1], title="v"; kwargs...)
    η_plt = contour(x, y[end:-1:1], η[:, end:-1:1], title="η"; kwargs...)

    plt1 = plot(u_plt, v_plt, η_plt, layout = (1,3))
    savefig(plt1, file)

end


function plot_2D_vorticity(phys, k_index, k, σmodes, geometry, mode_number, file)

    if geometry == Cartesian()
        y = phys[geometry].grid.y
        ylabel = "y (m)"
        Dy = phys[geometry].grid.D
    elseif geometry == Spherical()
        y = rad2deg.(grid.ϕ)
        ylabel = "ϕ (°)"
        Dy = phys[geometry].grid.D/phys[geometry].grid.a
    end 

    N = phys[geometry].grid.N
    Lx = 2*pi/k;
    Nx = N;
    x = LinRange(0, Lx, Nx+1)
  X, Y = meshgrid(x,y)
  
  uvec = σmodes[(k,geometry)][    1:N+1,   mode_number];
  vvec = σmodes[(k,geometry)][  N+2:2*N+2, mode_number];

    ζvec = im * k * vvec - Dy * uvec    
    ζ    = repeat(real(ζvec),1,Nx+1).*cos.(k*X) - repeat(imag(ζvec),1,Nx+1).*sin.(k*X);

    kwargs = (
        xlabel = "x",
        ylabel = ylabel,
          fill = true,
        levels = 20,
     linewidth = 0,
         color = :balance,
      colorbar = true,
          xlim = (    0, Lx),
         ylim = (y[end], y[1])
    )

    plt = contour(x, y[end:-1:1], ζ[:, end:-1:1], title="ζ"; kwargs...)

    savefig(plt, file)

end
