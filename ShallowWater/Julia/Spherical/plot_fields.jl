function plot_basic_state(y, U, E, file)

    p1 = plot( 
        y, 
        U,
        linewidth = 3,
        label = "U",
        title = "Velocity",
        xlims = (y[end], y[1])
        )
        
    p2 = plot(
        y, 
        E,
        linewidth = 3,
        label = "η",
        xlabel = "y",
        title = "Free-Surface",
        xlims = (y[end], y[1])
        )

    plt1 = plot(p1, p2, layout = (2,1))    
    savefig(plt1, file)

end

function plot_growth_rates(ks, σ, Nmodes, file)

    titles = ["1st mode", "2nd mode"]
    plt2 = plot()
    for cnt in 1:Nmodes
        plot!(plt2, 
        ks, 
        σ[cnt,:], 
        linewidth = 3,
        label=titles[cnt],
        title="Growth Rates",
        xlabel="k",
        xlims=(ks[1], ks[end])
        )
    end

    savefig(plt2,file)

end

function plot_1D_streamfunction(k_index, k, y, σmodes, mode_number, file)

      Ny = length(y) - 1
      Lx = 2*pi/k;
      Nx = Ny;
      dx = Lx/Nx;
       x = collect(0:dx:Lx);
    X, Y = meshgrid(x,y)
    
    # Divide eigenvector matrix into specified spaces: u,v,eta.
    uvec = σmodes[1:Ny+1,        mode_number, k_index];
    vvec = σmodes[Ny+2:2*Ny+2,   mode_number, k_index];
    ηvec = σmodes[2*Ny+1:3*Ny+1, mode_number, k_index];
    
    p1 = plot(
        y, 
        real(uvec),
        linewidth = 3, 
        label = "real",
        xlabel = "ϕ",
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

function plot_2D_streamfunction(k_index, k, y, σmodes, mode_number, file)

    Ny = length(y) - 1
    Lx = 2*pi/k;
    Nx = Ny;
    dx = Lx/Nx;
     x = collect(0:dx:Lx+dx/2);
  X, Y = meshgrid(x,y[end:-1:1])

    uvec = σmodes[1:Ny+1,        mode_number, k_index];
    vvec = σmodes[Ny+2:2*Ny+2,   mode_number, k_index];
    ηvec = σmodes[2*Ny+1:3*Ny+1, mode_number, k_index];

    u =     repeat(real(uvec),1,Nx+1).*cos.(k*X) -    repeat(imag(uvec),1,Nx+1).*sin.(k*X);
    v = -k.*repeat(imag(vvec),1,Nx+1).*cos.(k*X) - k.*repeat(real(vvec),1,Nx+1).*sin.(k*X);
    η =     repeat(real(ηvec),1,Nx+1).*cos.(k*X) -    repeat(imag(ηvec),1,Nx+1).*sin.(k*X);

    kwargs = (
        xlabel = "x",
        ylabel = "y",
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


function plot_2D_vorticity(Dy, k_index, k, y, σmodes, mode_number, file)

    Ly = params.Lϕ
    Ny = length(y) - 1
    Lx = 2*pi/k;
    Nx = Ny;
    dx = Lx/Nx;
     x = collect(0:dx:Lx+dx/2);
  X, Y = meshgrid(x,y[end:-1:1])

    uvec = σmodes[1:Ny+1,        mode_number, k_index];
    vvec = σmodes[Ny+2:2*Ny+2,   mode_number, k_index];

    ζvec = im * k * vvec - Dy * uvec    
    ζ    = repeat(real(ζvec),1,Nx+1).*cos.(k*X) - repeat(imag(ζvec),1,Nx+1).*sin.(k*X);

    kwargs = (
        xlabel = "x",
        ylabel = "y",
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
