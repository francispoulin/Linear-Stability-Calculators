using Plots
using IJulia

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

