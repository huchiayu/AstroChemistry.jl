function grain_recomb_H(Temp, psi)
    @assert (Temp > 0.0 && psi > 0.0)
    return 12.25e-14 / (1.0 + 8.074e-6 * psi^1.378 * (1.0 + 508.7*Temp^0.01586*psi^(-0.4723-1.102e-5*log(Temp))))
end

function grain_recomb_He(Temp, psi)
    return 5.572e-14 / (1.0 + 3.185e-7 * psi^1.512 * (1.0 + 5115*Temp^3.903e-7*psi^(-0.4956-5.494e-7*log(Temp))))
end

function grain_recomb_C(Temp, psi)
    return 45.58e-14 / (1.0 + 6.089e-3 * psi^1.128 * (1.0 + 433.1*Temp^0.04845*psi^(-0.8120-1.333e-4*log(Temp))))
end

function grain_recomb_Si(Temp, psi)
    return 2.166e-14 / (1.0 + 5.678e-8 * psi^1.874 * (1.0 + 43750*Temp^1.635e-6*psi^(-0.8964-7.538e-5*log(Temp))))
end
