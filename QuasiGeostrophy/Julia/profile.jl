function profile(y, params);

    Ψ  = -params.Uj * tanh.(y / params.Lj);
    U  =  params.Uj * sech.(y / params.Lj).^2;

    return U,Ψ
end