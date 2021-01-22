function profile(y, params);

    Δη = params.Uj * params.Lj * params.f₀ / params.g

    U  =  params.Uj * sech.(y / params.Lj).^2;
    E  = -       Δη * tanh.(y / params.Lj);

    return U,E
end