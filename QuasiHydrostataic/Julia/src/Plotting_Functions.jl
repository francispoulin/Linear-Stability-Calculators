using Plots
using ColorSchemes
using Printf

pyplot()

function meshgrid(x, y)
    x = reshape(x, 1, length(x))
    y = reshape(y, length(y), 1)
    X = repeat(x, length(y), 1)
    Y = repeat(y, 1, length(x))
    return X, Y
end

function plot_growth_slice(filenc, filejson, fileplot)

    # read json file
    json_string = JSON3.read(read(filejson, String))

    Umax = json_string.jet.Umax
    L    = json_string.jet.L
    fz   = json_string.physics.fz

    # read NetCDF file
    ωs, modes, ks, ms, iEigs, ys = load_spectrum(filenc)

    print("\n--> Plotting the growth rates for a fixed k versus m in ", fileplot, "\n")

    growth = zeros(size(ωs))
    for (ik, k) in enumerate(ks)
        for (im, m) in enumerate(ms)
            for iEig in iEigs
                growth[ik, im, iEig] = max(0, imag(ωs/fz)[ik, im, iEig])
            end
        end
    end

    # scaling wavenumber for plotting
    m = L/(2*π)*ms
    plt = plot()
    for iEig in iEigs
        plot!(plt, m, growth[1,:,iEig],
                label = @sprintf("%d", iEig),
            linewidth = 2,
            xlim = (m[1], m[end]),
               xlabel = "L m/(2*π)",
                title = "Growth rates for k = 0",
               legend = :topright)
    end
    savefig(plt, fileplot)
end 

function plot_growth(filenc, filejson, fileplot)

    # read json file
    json_string = JSON3.read(read(filejson, String))

    Umax = json_string.jet.Umax
    L    = json_string.jet.L
    fz   = json_string.physics.fz

    # read NetCDF file
    ωs, modes, ks, ms, iEigs, ys = load_spectrum(filenc)

    print("\n--> Plotting the growth rates versus (k, m) in ", fileplot, "\n")

    growth = zeros(size(ωs)[1:2])
    for (ik, k) in enumerate(ks)
        for (im, m) in enumerate(ms)
            growth[ik, im] = max(0, imag(ωs/fz)[ik, im, 1])
        end
    end

    # scaling wavenumber for plotting
    k = L/(2*π)*ks
    m = L/(2*π)*ms

    kwargs = (
        xlabel = "k",
        ylabel = "m",
          fill = true,
        levels = 20,
     linewidth = 0,
         color = :haline, 
      colorbar = true,
          xlim = (k[1], k[end]),
          ylim = (m[1], m[end])
    )

    plt = contour(k, m, growth[:,:,1]', title="Growth rates for 1st mode"; kwargs...)
    savefig(plt, fileplot)
end 

function plot_modes_1D(filenc, filejson, fileplot, Neigs)

        # read json file
        json_string = JSON3.read(read(filejson, String))

        Umax = json_string.jet.Umax
        L    = json_string.jet.L
        fz   = json_string.physics.fz
        Ny   = json_string.grid.Ny
        Ly   = json_string.grid.Ly
        
        # read NetCDF file
        ωs, modes, ks, ms, iEigs, ys = load_spectrum(filenc)
    
        print("\n--> Plotting 1D modes in the file ", fileplot, "\n")

        # scaling wavenumber for plotting
        y = ys[1:Ny+1]
        k = L/(2*π)*ks
        m = L/(2*π)*ms

        ik = 1    # most unstable mode
        
        for iEig in collect(1:Neigs)
            tmp              = findmax(imag(ωs/fz)[ik, :, iEig])
            max_growth, indm = tmp[1], tmp[2]

            @printf("\n iEig = %2d  k = %10.6f  m = %10.6f  growth = %10.6f", 
                    iEig, float(k[ik]), float(m[indm]), max_growth)

            vvec          = zeros(ComplexF64, Ny+1)
            uvec          = modes[ik, indm, 1:Ny+1, iEig] 
            vvec[2:end-1] = modes[ik, indm, Ny+2:2*Ny, iEig]
            bvec          = modes[ik, indm, 2*Ny+1:3*Ny+1, iEig] 

            u_plt = plot(y/ L, real(uvec), line=:solid, color=:blue, linewidth=2, label="Re(u)")
            plot!(u_plt, y/ L, imag(uvec), line=:solid, color=:red,  linewidth=2, label="Im(u)")
            v_plt = plot(y/ L, real(vvec), line=:solid, color=:blue, linewidth=2, label="Re(v)")
            plot!(v_plt, y/ L, imag(uvec), line=:solid, color=:red,  linewidth=2, label="Im(v)")
            b_plt = plot(y/ L, real(bvec), line=:solid, color=:blue, linewidth=2, label="Re(b)")
            plot!(b_plt, y/ L, imag(bvec), line=:solid, color=:red,  linewidth=2, label="Im(b)")
            plt = plot(u_plt, v_plt, b_plt, layout=(3,1), size=(800, 600))  
            file = string("modes_1D_iEig", iEig, ".png")
            savefig(plt, file)
        end
end

function plot_modes_2D(filenc, filejson, fileplot, Neigs)

    # read json file
    json_string = JSON3.read(read(filejson, String))

    Umax = json_string.jet.Umax
    L    = json_string.jet.L
    fz   = json_string.physics.fz
    Ny   = json_string.grid.Ny
    Ly   = json_string.grid.Ly
    
    # read NetCDF file
    ωs, modes, ks, ms, iEigs, ys = load_spectrum(filenc)

    print("\n--> Plotting 2D modes in the file ", fileplot, "\n")

    # scaling wavenumber for plotting
    y = ys[1:Ny+1]
    k = L/(2*π)*ks
    m = L/(2*π)*ms

    ik = 1    # most unstable mode
    iEig = 1
    indm = findmax(imag(ωs/fz)[ik, :, iEig])[2]
    Nz = Ny
    Lz = 2*π/ms[indm]
    dz = Lz/Nz
    z = collect(0:dz:Lz)
    Y,Z = meshgrid(y, z)

    kwargs = (
        xlabel = "y",
        ylabel = "z",
          fill = true,
        levels = 20,
     linewidth = 0,
         color = :balance,
      colorbar = true,
#          xlim = (k[1], k[end]),
#          ylim = (m[1], m[end])
    )

    for iEig in collect(1:Neigs)
        tmp              = findmax(imag(ωs/fz)[ik, :, iEig])
        max_growth, indm = tmp[1], tmp[2]

        @printf("\n iEig = %2d  k = %10.6f  m = %10.6f  growth = %10.6f", 
                iEig, float(k[ik]), float(m[indm]), max_growth)

        # 1D structures
        uvec            = zeros(ComplexF64, 1, Ny+1)
        vvec            = zeros(ComplexF64, 1, Ny+1)
        bvec            = zeros(ComplexF64, 1, Ny+1)
        uvec[1,:]       = modes[ik, indm, 1:Ny+1, iEig] 
        vvec[1,2:end-1] = modes[ik, indm, Ny+2:2*Ny, iEig]
        bvec[1,:]       = modes[ik, indm, 2*Ny+1:3*Ny+1, iEig] 

        # 2D structures
        U = real(repeat(uvec, Nz+1, 1) .* exp.(1im * ms[indm] * Z))
        V = real(repeat(vvec, Nz+1, 1) .* exp.(1im * ms[indm] * Z))
        B = real(repeat(bvec, Nz+1, 1) .* exp.(1im * ms[indm] * Z))

        u_plt = contour(Y/1e3, Z/1e3, U', title="u"; kwargs...)
        v_plt = contour(Y/1e3, Z/1e3, V', title="v"; kwargs...)
        b_plt = contour(Y/1e3, Z/1e3, B', title="b"; kwargs...)
        plt = plot(u_plt, v_plt, b_plt, layout=(1,3), size=(1600, 400))  
        file = string("modes_2D_iEig", iEig, ".png")
        savefig(plt, file)
    end
end

#=
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
=#

#string.(fieldnames(typeof(gridC)))

# Plot vertical slices as well?E