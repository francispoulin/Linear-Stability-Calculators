using LinearAlgebra

"""
Inertial instability version of build_matrix
"""

function build_matrix(k, m, grid, physics, jet)

    N2     =  physics.N^2        
    Z      = zeros(grid.Ny+1, grid.Ny+1)
    I      =  diagm(0 => ones(grid.Ny+1))
    Q      =  diagm(0 =>   physics.fz .- jet.dU)
    By     =  diagm(0 => - physics.fy .* jet.dU)
    U      = k * diagm(0=>jet.U) - 1im * physics.ν * m^2 * I
    fyomDy = physics.fy / m * grid.D

    # Set up eigenvalue problem: [u, v, b]
    A = [                             U                (1im*Q + fyomDy)[:, 2:end-1]   -1im*k/m*I; 
        -(1im*physics.fz*I + fyomDy)[2:end-1,:]                  U[2:end-1,2:end-1]  -1/m*grid.D[2:end-1,:] ;
        1im*N2*k/m*I                             (N2/m*grid.D - 1im* By)[:,2:end-1]            U];

    return A
end

function compute_spectrum(A, modes, Neigs)

    λ_in, λvec_in = eigen(A);
            index = sortperm(imag(λ_in), rev=true);
             λvec = λvec_in[:,index];
                λ = λ_in[index];
     
               ωs = λ[1:Neigs];
            modes = λvec[:,1:Neigs];
    
        return ωs, modes
    end

    # Must modify this for this particular problem