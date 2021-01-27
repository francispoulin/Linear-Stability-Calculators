using LinearAlgebra: eigen, sortperm

function compute_spectrum(A, B, k, params, Nmodes)

λ_in, λvec_in = eigen(A,B);
        index = sortperm(imag(λ_in), rev=true);
         λvec = λvec_in[:,index];
            λ = λ_in[index];
 
       growth = imag(λ[1:Nmodes]) * k;
         freq = real(λ[1:Nmodes]) * k;
        modes = λvec[:,1:Nmodes];

    return growth, freq, modes
end