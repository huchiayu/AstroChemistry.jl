@with_kw struct Par{NPIX,N_neq,T}
    nH::T
    temp::T
    #xelec::T
    ξ::T
    IUV::T
    Zp::T
    NH ::SVector{NPIX,T} = zeros(NPIX)
    NH2::SVector{NPIX,T} = zeros(NPIX)
    NCO::SVector{NPIX,T} = zeros(NPIX)
    NC ::SVector{NPIX,T} = zeros(NPIX)
    xneq::SVector{N_neq,T} = []
end


function calc_coeff!(coeff, par::Par, xelec, net::Network{N_spec, N_reac, N_integrate, N_neq, T}) where{N_spec, N_reac, N_integrate, N_neq, T}
    @unpack nH, temp, NH, NH2, NCO, NC, ξ, IUV, Zp = par
    @unpack typ, alpha, beta, gamma, grRec, SupTh, ssC, iCion = net
    Av = NH .* (5.35e-22 * Zp)
    ξp = ξ / ξ0;
    for i in 1:N_reac
        if typ[i] == 1
            coeff[i] = IUV * alpha[i] * mean(@. exp(-gamma[i] * Av))
        elseif typ[i] == 2
            coeff[i] = alpha[i] * ξp
        elseif typ[i] == 3
            coeff[i] = alpha[i] * ξp * (temp / 300.)^beta[i] * gamma[i] / (1.0 - omega)
        elseif typ[i] == 0
            coeff[i] = alpha[i] * (temp / 300.)^beta[i] * exp(-gamma[i] / temp)
        end
        #overwrite
        if SupTh == true && typ[i] == 5 && NH2 < 4e0
            #=
            m_reduce = 1.714  #m_C+ = 12, m_H2 = 2, 12*2/(12+2) = 1.714
            vA = 3.3 #km/sec
            temp_eff = temp + 40.3825 * m_reduce * vA^2
            coeff[i] = alpha[i] * (temp_eff / 300.)^beta[i] * exp(-gamma[i] / temp_eff)
            if isnan(coeff[i])
                @show "nan!!!", i, coeff[i], temp_eff
            end
            =#
            coeff[i] = alpha[i] * (800. / 300.)^beta[i] * exp(-gamma[i] / 800.)
        end
        if isnan(coeff[i])
            @show i
        end
    end
    #=
    if temp <= 261.
        coeff[re[5643]] = 3.5e-11
    else
        coeff[re[5643]] = 1.77e-11 * exp(178.0/temp)
    end
    lnT = log(temp)
    coeff[re[6162]] = 2.753e-14 * (315614.0/temp)^1.5 * (1.0 + (115188.0/temp)^0.407)^(-2.242)
    coeff[re[6166]] = 1e-11 * temp^(-0.5) * (11.19 - 1.676*lnT - 0.2852*lnT^2 + 0.04433*lnT^3)
    coeff[re[3431]] = 1.4e-9 * (temp/300.)^(-0.5)
    coeff[re[1158]] = 7.0e-8 * (temp/300.)^(-0.5)
    coeff[re[1293]] = 1.08e-7 * (temp/300.)^(-0.5)
    coeff[re[1295]] = 6.02e-8 * (temp/300.)^(-0.5)
    coeff[re[1296]] = 2.58e-7 * (temp/300.)^(-0.5)
    =#
    #photodissociation
    coeff[end]   = IUV * kH2diss * mean(@. exp(-gamma_H2 * Av) * fH2selfshield(NH2))
    coeff[end-1] = IUV * kCOdiss * mean(@. exp(-gamma_CO * Av) * fCOselfshield(NCO,NH2))
    #coeff[end-1] = fCOselfshield(1.,1.)

    #grain processes
    #coeff[end-2] = kdust * (temp / 100.)^0.5 * Zp
    coeff[end-2] = get_kdust(temp) * Zp
    psi = 1e20 #large value will turn off grain recombination
    #@show grRec
    if grRec == true
        if xelec > 1e-20
            psi = 1.7 * (IUV+1e-20) * mean(@. exp(-1.87*Av)) * sqrt(temp) / (nH * xelec)
        end
        psi = psi < 1e20 ? psi : 1e20
    end

    if psi <= 0.0
        println("IUV=", IUV, "  AV=", Av, "  temp=", temp, "  nH=", nH, "  xelec=", xelec)
        error("psi<=0")
    end
    coeff[end-3] = grain_recomb_H( temp, psi) * Zp
    coeff[end-4] = grain_recomb_He(temp, psi) * Zp
    coeff[end-5] = grain_recomb_C( temp, psi) * Zp
    coeff[end-6] = grain_recomb_Si(temp, psi) * Zp

    if ssC == true
        coeff[iCion] *= mean(@. fCselfshield(NC,NH2))
    end
    #=
    if SupTh == true && NH2 < 4e20
        m_reduce = 1.714  #m_C+ = 12, m_H2 = 2, 12*2/(12+2) = 1.714
        vA = 3.3 #km/sec
        temp_eff = temp + 40.3825 * vA^2
        coeff[iSupTh] = alpha[iSupTh] * (temp_eff / 300.)^beta[iSupTh] * exp(-gamma[iSupTh] / temp_eff)
    end
    =#
    #if fCselfshield(NC,NH2) < 0.1 @show NC, NH2, NCO end
end

function get_kdust(temp)
    deff=1.0
    dust_to_gas_ratio = 1.0
    tdust = 15.
    tdust2 = tdust * 1e-2
    fa     = 1.0 / (1.0 + 1e4 * exp(-6e2 / tdust))

    temp2  = temp * 1e-2
    stick  = 1.0 / (1.0 + 0.4 * (temp2 + tdust2)^0.5 + 0.2 * temp2 + 0.08 * temp2^2)
    ch7    = dust_to_gas_ratio * deff * 3e-17 * temp2^0.5
    Rdust  = ch7 * fa * stick
end

function init_abund(abund, Zp, xneq, abtot, net)
    #abund[dict["C+"]] = abC_s * Zp
    #abund[dict["O"]]  = abO_s * Zp
    calc_abund_derived(abund, Zp, xneq, abtot, net)
end

function calc_abund_derived(abund, Zp, xneq, abtot::AbundTotal, net::Network{N_spec, N_reac, N_integrate, N_neq, T}) where{N_spec, N_reac, N_integrate, N_neq, T}

    @unpack abC_s, abO_s, abSi_s, abS_s, abN_s, abFe_s, abMg_s = abtot
    @unpack fac_H, fac_He, fac_C, fac_O, fac_Si, fac_S, fac_N, fac_Mg, fac_Fe, charge, ineq, iH, iHe, iC, iO, iSi, iS, iN, iMg, iFe, ielec = net

    if N_neq > 0
        for i in 1:N_neq
            abund[ineq[i]] = xneq[i]
        end
    end

    sH = sHe = sC = sO = sSi = sS = sN = sMg = sFe = selec = 0.0
    for i in 1:N_spec
        (iH  > 0) && (i != iH)  ? (sH  += abund[i] * fac_H[i])  : nothing
        (iHe > 0) && (i != iHe) ? (sHe += abund[i] * fac_He[i]) : nothing
        (iC  > 0) && (i != iC)  ? (sC  += abund[i] * fac_C[i])  : nothing
        (iO  > 0) && (i != iO)  ? (sO  += abund[i] * fac_O[i])  : nothing
        (iSi > 0) && (i != iSi) ? (sSi += abund[i] * fac_Si[i]) : nothing
        (iS  > 0) && (i != iS)  ? (sS  += abund[i] * fac_S[i])  : nothing
        (iN  > 0) && (i != iN)  ? (sN  += abund[i] * fac_N[i])  : nothing
        (iMg > 0) && (i != iMg) ? (sMg += abund[i] * fac_Mg[i]) : nothing
        (iFe > 0) && (i != iFe) ? (sFe += abund[i] * fac_Fe[i]) : nothing
        (ielec  > 0) && (i != ielec) ? (selec += abund[i] * charge[i]) : nothing
    end

    iH  > 0 ? abund[iH]  = 1.0 - sH : nothing
    iHe > 0 ? abund[iHe] = XHe - sHe : nothing
    iC  > 0 ? abund[iC]  = abC_s * Zp - sC : nothing
    iO  > 0 ? abund[iO]  = abO_s * Zp - sO : nothing
    iS  > 0 ? abund[iS]  = abS_s * Zp - sS : nothing
    iSi > 0 ? abund[iSi] = abSi_s * Zp - sSi : nothing
    iN  > 0 ? abund[iN]  = abN_s * Zp - sN : nothing
    iMg > 0 ? abund[iMg] = abMg_s * Zp - sMg : nothing
    iFe > 0 ? abund[iFe] = abFe_s * Zp - sFe : nothing
    ielec > 0 ? abund[ielec] = selec : nothing

end

function solve_equilibrium_abundances(abund, dtime, par::Par, abtot::AbundTotal, net::Network{N_spec, N_reac, N_integrate, N_neq, T}) where{N_spec, N_reac, N_integrate, N_neq, T}

    @unpack typ, ir1, ir2, ip1, ip2, ip3, ip4, pro3, pro4, alpha, beta, gamma, idx_integrate, ielec, iH2 = net

    abund_dot = zeros(N_spec)
    abund_final = zeros(N_spec)
    abund_con = zeros(N_integrate)

    coeff = zeros(N_reac)
    reac_rates = zeros(N_reac)

    @unpack nH, Zp, xneq = par

    max_abund_inv = zeros(N_spec)
    get_max_abundance!(max_abund_inv, Zp, abtot, net)

    #we use closure s.t. large arrays like coeff & reac_rates can be updated in-place
    function calc_abund_dot_closure(abund_dot_int, abund_int, ppp, t)

        abund_dot .= 0.0 #in-place for better performance

        abund[idx_integrate] .= abund_int
        calc_abund_derived(abund, Zp, xneq, abtot, net)

        calc_coeff!(coeff, par, abund[ielec], net) #zero allocation

        for i in 1:N_reac
            if typ[i] == 0
                reac_rates[i] = coeff[i] * abund[ir1[i]] * abund[ir2[i]] * nH  #2-body process
            elseif typ[i] == 3
                reac_rates[i] = coeff[i] * abund[ir1[i]] * abund[iH2]  #CR-induced photon process
            elseif typ[i] == 4
                reac_rates[i] = coeff[i] * abund[ir1[i]] * nH  #grain processes
            else
                #typ[i] == 1 or 2
                reac_rates[i] = coeff[i] * abund[ir1[i]]  #CR or photo process
            end

            #destruction
            if ir1[i] != 0
                abund_dot[ir1[i]] -= reac_rates[i]
            end
            if ir2[i] != 0
                abund_dot[ir2[i]] -= reac_rates[i]
            end

            #formation
            if ip1[i] != 0
                abund_dot[ip1[i]] += reac_rates[i]
            end
            if ip2[i] != 0
                abund_dot[ip2[i]] += reac_rates[i]
            end
            if pro3 == true
                if ip3[i] != 0
                    abund_dot[ip3[i]] += reac_rates[i]
                end
            end
            if pro4 == true
                if ip4[i] != 0
                    abund_dot[ip4[i]] += reac_rates[i]
                end
            end
        end #zero allocation
        abund_dot_int .= abund_dot[idx_integrate]
    end
    rtol, atol, ofbtol = 1e-4, 1e-7, 1e-2
    #tend = SEC_PER_YEAR * 1e10 / Zp / nH
    tend = SEC_PER_YEAR * dtime
    tspan = (0.0, tend)
    abund_con .= abund[idx_integrate] #can't use idx_integrate
    prob = ODEProblem(calc_abund_dot_closure, abund_con, tspan)
    sol = solve(prob,
        alg_hints=[:stiff],
        reltol=rtol, abstol=atol,
        #isoutofdomain=(y,p,t)->any(x->(x<0.0-ofbtol||x>1.0+ofbtol),abund.*max_abund_inv),
        isoutofdomain=(y,p,t)->speciesoutofbound(abund, max_abund_inv, ofbtol),
        maxiters=1e3,
        save_everystep=false);
    abund_final[idx_integrate] .= sol.u[end]
    calc_abund_derived(abund_final, Zp, xneq, abtot, net)

    abund .= abund_final
    speciesoutofbound(abund, max_abund_inv, Inf) #fixing OFB error
    return sol.retcode, reac_rates
end


function speciesoutofbound(abund, max_abund_inv, ofbtol)

    #for i in 1:N_spec
    for i in eachindex(abund)
        x = abund[i] * max_abund_inv[i]
        if x < 0.0 - ofbtol
            #println("species out of bound, species=", all_species[i], "  abund=", abund[i])
            return true
        elseif x < 0.0
            #println("fixing small error, species=", all_species[i], "  abund=", abund[i])
            abund[i] = 0.0
        end
        if x > 1.0 + ofbtol
            #println("species out of bound, species=", all_species[i], "  abund=", abund[i])
            return true
        elseif x > 1.0
            #println("fixing small error, species=", all_species[i], "  abund=", abund[i])
            abund[i] = 1.0 / max_abund_inv[i]
        end
    end
    return false
end

function get_max_abundance!(max_abund_inv, Zp, abtot::AbundTotal, net::Network)
    #update in-place
    @unpack abC_s, abO_s, abSi_s, abS_s, abN_s, abFe_s, abMg_s = abtot
    @unpack fac_H, fac_He, fac_C, fac_O, fac_S, fac_Si, fac_N, fac_Mg, fac_Fe = net
    for i in eachindex(max_abund_inv)
        val = 0
        val = max(fac_H[i], fac_He[i]/XHe)
        val = abC_s  > 0 ? max(fac_C[i] /(abC_s*Zp ), val) : val
        val = abO_s  > 0 ? max(fac_O[i] /(abO_s*Zp ), val) : val
        val = abS_s  > 0 ? max(fac_S[i] /(abS_s*Zp ), val) : val
        val = abSi_s > 0 ? max(fac_Si[i]/(abSi_s*Zp), val) : val
        val = abN_s  > 0 ? max(fac_N[i] /(abN_s*Zp ), val) : val
        val = abMg_s > 0 ? max(fac_Mg[i]/(abMg_s*Zp), val) : val
        val = abFe_s > 0 ? max(fac_Fe[i]/(abFe_s*Zp), val) : val
        max_abund_inv[i] = val
    end
end
