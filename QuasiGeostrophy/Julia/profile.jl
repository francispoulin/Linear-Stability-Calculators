function profile(y, params);

    Uj = params.Uj
    Lj = params.Lj
    Fr = params.Fr

    Ψ  = -Uj * tanh.(y/Lj);
    U  =  Uj * sech.(y/Lj).^2;
    Qy = Fr*U .- Uj*(4 * sinh.(y/Lj).^2 .- 2) .* sech.(y/Lj).^4;

    return Qy,U,Ψ
end

# Q    = Ψ_yy - Fr Ψ
# Qy   = Ψ_yyy - Fr Ψ_y = - U_yy + Fr U = Fr U - U_yy
# U_yy = Uj * ( 4*tanh.(y/Lj).^2*sech.(y/Lj).^2 - sech.(y/Lj).^4) 