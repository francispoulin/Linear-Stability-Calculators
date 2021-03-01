using Plots
using IJulia

# choose y label to be either y or ϕ, depending on the geometry
function plot_basic_state(y, solution, file)

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

