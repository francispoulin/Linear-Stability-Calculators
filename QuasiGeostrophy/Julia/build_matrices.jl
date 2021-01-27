using LinearAlgebra

function build_AB(k, dQdy, U, Ψ, Dy, Dy2, y, params)

    k2 = k*k;
    Fr = params.Fr;
    Ny = length(U)-1;

    B = Dy2[2:end-1,2:end-1] - (k2 + Fr)*I;
    A = diagm(U[2:end-1]) * B + diagm(dQdy[2:end-1]);

    return A, B
end