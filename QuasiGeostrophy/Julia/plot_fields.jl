function plot_basic_state(y, dQdy, U, Ψ, Ly, file)

    p1 = plot(
        y, 
        Ψ,
        linewidth = 3,
        label = "Ψ",
        xlabel = "y",
        title = "Streamfunction",
        xlims = (-Ly/2, Ly/2)
        )

    p2 = plot(
        y, 
        U,
        linewidth = 3,
        label = "U",
        title = "Velocity",
        xlims = (-Ly/2, Ly/2)
        )
            
    p3 = plot(
        y, 
        dQdy,
        linewidth = 3,
        label = "dQ/dy",
        xlabel = "y",
        title = "dQ/dy",
        xlims = (-Ly/2, Ly/2)
        )
    
    plt1 = plot(p1, p2, p3, layout = (3,1))    
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
        title="Growth Rates of QG Model",
        xlabel="k",
        xlims=(ks[1], ks[end]),
        xticks=(0:0.25:ks[end]),
        gridlinewidth=2.0
        )
    end

    savefig(plt2,file)

end

function plot_1D_fields(k_index, k, Ny, y, σmodes, mode_number, file)

      Ly = params.Ly
      Lx = 2*pi/k;
      Nx = Ny;
      dx = Lx/Nx;
       x = collect(0:dx:Lx);
    X, Y = meshgrid(x,y)
    
    # Divide eigenvector matrix into specified spaces: u,v,eta.
    ψvec = σmodes[:, mode_number, k_index];
    
    plt = plot(
        y, 
        real(ψvec),
        linewidth = 3, 
        label = "real",
        xlabel = "y",
        title = "Streamfunction",
        xlims = (-Ly/2, Ly/2)
        )
    plot!(
        plt,
        y,
        linewidth = 3,
        imag(ψvec),
        label = "imag")
    savefig(plt, file)
    #display(plt)
    
end

function plot_2D_fields(k_index, k, Ny, y, σmodes, mode_number, file)

    Ly = params.Ly
    Lx = 2*pi/k;
    Nx = Ny;
    dx = Lx/Nx;
     x = collect(0:dx:Lx+dx/2);
  X, Y = meshgrid(x,y[end:-1:1])

    ψvec = σmodes[1:Ny+1,        mode_number, k_index];

    ψ = repeat(real(ψvec),1,Nx+1).*cos.(k*X) - repeat(imag(ψvec),1,Nx+1).*sin.(k*X);

    kwargs = (
        xlabel = "x",
        ylabel = "y",
          fill = true,
        levels = 20,
     linewidth = 0,
         color = :balance,
      colorbar = true,
          xlim = (    0, Lx),
         ylim =  (-Ly/2, Ly/2)
    )

    plt = contour(x, y[end:-1:1], ψ[:, end:-1:1], title="ψ"; kwargs...)

    savefig(plt, file)

end