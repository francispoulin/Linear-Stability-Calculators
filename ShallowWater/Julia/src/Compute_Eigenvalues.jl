using LinearAlgebra

mutable struct Build_Matrix{F, S, P}
           k::F
    solution::S
        phys::P
end

"""
Cartesian version
"""
function Build_Matrix(geometry::Cartesian;
    k,
    solution,
    phys
    )

    A = zeros(phys.grid.N)

    return Build_Matrix(A)
end

"""
Spherical version
"""
function Build_Matrix(geometry::Spherical;
    k,
    solution,
    phys
    )

    A = zeros(phys.grid.N)

    return Build_Matrix(A)
end

