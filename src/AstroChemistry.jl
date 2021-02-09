#push!(LOAD_PATH, pwd())
module AstroChemistry

species_noneq = []
const N_neq = length(species_noneq)
const NONEQ = N_neq > 0
const grRec = true

export solve_equilibrium_abundances, calc_abund_derived, init_abund
export calc_coeff!
export kH2diss, kdust
export dict
export iH2,iCO,iC
export fac_H, fac_C, fac_O, charge
export abC_s, abO_s
export N_spec, N_reac
export Par
export fH2selfshield, fCOselfshield

using DelimitedFiles
using Printf
using StaticArrays
using DifferentialEquations
#using OrdinaryDiffEq
#using NLsolve
#using PyPlot  #can't use Interact with it
#using Plots
#plotlyjs()
#using Interpolations
using Parameters
using Statistics

#using Sundials #for CVODE_BDF

include("OctTree.jl")
include("GadgetReader.jl")


include("shielding_functions.jl")
include("init_umist.jl")
include("grain_recomb.jl")


all_species =
[
"H", "H-", "H2", "H+", "H2+", "H3+", "e-", "He", "He+", "HeH+",
"C", "C+", "CO", "HCO+",
"O", "O+", "OH", "OH+", "H2O+", "H3O+", "H2O",
"O2", "CO+", "O2+",
"CH2", "CH2+", "CH", "CH+", "CH3+",
"Si+", "Si"
#"SiO", "SiO+", "SiOH+",
#"S", "S+",
]

const N_spec = size(all_species,1)

all_species_2 = vcat(all_species , ["", "PHOTON", "CRP", "CRPHOT"]);

reac_num=[]
for i in 1:size(umist,1)
    #println(i)
    n, t, R, P, a, b, c = read_umist_line(umist[i,:])
    if any(R[1].==all_species_2) && any(R[2].==all_species_2) && any(P[1].==all_species_2) && any(P[2].==all_species_2) && any(P[3].==all_species_2) && any(P[4].==all_species_2)
        #print_reaction(n, t, R, P, a, b, c)
        push!(reac_num,n)
    end
end

reduced_network = umist[reac_num,:]

#add 7 more reactions (H2 & CO dissociation & grain-processes)
more_reactions = readdlm("src/user-rates.txt",':')
reduced_network = vcat(reduced_network, more_reactions);
const N_reac = size(reduced_network, 1)
#@show N_reac, N_spec

dict = Dict(all_species .=> collect(1:N_spec) );
dict[""] = dict["PHOTON"] = dict["CRP"] = dict["CRPHOT"] = 0

#this need to be after dict is constructed
include("constants.jl")


charge_a = @. Int32(getindex(breakdown_species(all_species) , 1))
elements = @. getindex(breakdown_species(all_species) , 2)
counts = @. getindex(breakdown_species(all_species) , 3)

function get_unique_elements(elements)
    unique_elements = String[]
    for i in 1:N_spec
        for j in eachindex(elements[i])
            push!(unique_elements, elements[i][j])
        end
    end
    return unique!(unique_elements)
end

function get_idx_integrate(elements)
    unique_elements = get_unique_elements(elements)
    #const N_uni_elem = length(unique_elements_a)
    #const uni_elem = SVector{N_uni_elem}(unique_elements_a)
    #N_conserve = length(unique_elements) + 1  #+1 for electron

    #species_derived = ["e-", "H", "He", "C", "O", "Si", "S"]
    species_derived = unique_elements
    push!(species_derived, "e-")
    #@show species_derived

    if NONEQ
        for i in 1:N_neq
            push!(species_derived, species_noneq[i])
        end
    end
    #@show species_derived
    #@show ineq
    #@show all_species
    idx_derived = sort([dict[species_derived[i]] for i in eachindex(species_derived)])

    idx_integrate_a = filter!(s -> all(s .!= idx_derived), collect(1:N_spec))
    return idx_integrate_a, species_derived
end

idx_integrate_a, species_derived = get_idx_integrate(elements)
const N_integrate = length(idx_integrate_a)
const idx_integrate = SVector{N_integrate}( idx_integrate_a )

if NONEQ
    const ineq = SVector{N_neq}([dict[species_noneq[i]] for i in 1:N_neq])
end

function construct_index_arrays(reduced_network)

    num_a = zeros(Int32, N_reac)
    typ_a = zeros(Int32, N_reac)

    ir1_a = zeros(Int32, N_reac)
    ir2_a = zeros(Int32, N_reac)
    ip1_a = zeros(Int32, N_reac)
    ip2_a = zeros(Int32, N_reac)
    ip3_a = zeros(Int32, N_reac)
    ip4_a = zeros(Int32, N_reac)

    alpha_a = zeros(N_reac)
    beta_a  = zeros(N_reac)
    gamma_a = zeros(N_reac)

    for i in 1:N_reac
        num_a[i], t, R, P, alpha_a[i], beta_a[i], gamma_a[i] = read_umist_line(reduced_network[i,:])
        if any(R .== "C+") && any(R .== "e-")
            #alpha_a[i] = 4.67e-12
            #gamma_a[i] = 3.76
            #@show alpha_a[i], beta_a[i], gamma_a[i]
        end
        if any(R .== "H-") && any(R .== "PHOTON")
            #alpha_a[i] = 55.8e-10
            #gamma_a[i] = 3.76
            #@show alpha_a[i], beta_a[i], gamma_a[i]
        end
        if any(R .== "PHOTON") && any(R .== "C") #C photoionization
            iCion = i
            #@show alpha_a[i], gamma_a[i]
        end

        if any(R .== "PHOTON")
            typ_a[i] = 1
        elseif any(R .== "CRP")
            typ_a[i] = 2
        elseif any(R .== "CRPHOT")
            typ_a[i] = 3
        else
            typ_a[i] = 0
        end

        #overwrite if it's a grain-process
        if t == "GR"
            typ_a[i] = 4
        end

        #overwrite if it's an ion-neutral process
        if SupTh == true && t == "IN"
            if any(R .== "C+") && any(R .== "H2") && any(P .== "CH+") && any(P .== "H")
                typ_a[i] = 5
                @show i,R,P
            end
        end

        ir1_a[i], ir2_a[i] = dict[R[1]], dict[R[2]]
        ip1_a[i], ip2_a[i], ip3_a[i], ip4_a[i] = dict[P[1]], dict[P[2]], dict[P[3]], dict[P[4]]
    end
    return typ_a, ir1_a, ir2_a, ip1_a, ip2_a, ip3_a, ip4_a, alpha_a, beta_a, gamma_a
end
typ_a, ir1_a, ir2_a, ip1_a, ip2_a, ip3_a, ip4_a, alpha_a, beta_a, gamma_a = construct_index_arrays(reduced_network)

const typ = SVector{N_reac}(typ_a); const ir1 = SVector{N_reac}(ir1_a); const ir2 = SVector{N_reac}(ir2_a);
const ip1 = SVector{N_reac}(ip1_a); const ip2 = SVector{N_reac}(ip2_a); const ip3 = SVector{N_reac}(ip3_a);
const ip4 = SVector{N_reac}(ip4_a); const pro4 = any(ip4_a .!= 0); const pro3 = any(ip3_a .!= 0);

const alpha = SVector{N_reac}(alpha_a);
const beta  = SVector{N_reac}(beta_a);
const gamma = SVector{N_reac}(gamma_a);

function setup_conservation_arrays(elem)
    idx_spec = Int32[]
    count_spec = Int32[]
    fac = zeros(Int32, N_spec)

    for i in eachindex(elements)
        for j in eachindex(elements[i])
            #@show i, j, elements[i][j]
            if elements[i][j] == elem
                push!(idx_spec, i)
                push!(count_spec, counts[i][j])
                fac[i] = counts[i][j]
                #@show counts[i][j], all_species[i]
            end
        end
    end
    fac
end


const fac_H  = SVector{N_spec}(setup_conservation_arrays("H"))
const fac_He = SVector{N_spec}(setup_conservation_arrays("He"))
const fac_C  = SVector{N_spec}(setup_conservation_arrays("C"))
const fac_O  = SVector{N_spec}(setup_conservation_arrays("O"))
const fac_Si = SVector{N_spec}(setup_conservation_arrays("Si"))
const fac_S  = SVector{N_spec}(setup_conservation_arrays("S"))
const fac_N  = SVector{N_spec}(setup_conservation_arrays("N"))
const fac_Mg = SVector{N_spec}(setup_conservation_arrays("Mg"))
const fac_Fe = SVector{N_spec}(setup_conservation_arrays("Fe"))
const charge = SVector{N_spec}(charge_a);


@with_kw struct Par{NPIX,T}
    nH::T
    temp::T
    #xelec::T
    ξ::T
    IUV::T
    Zp::T
    NH ::SVector{NPIX,T}
    NH2::SVector{NPIX,T}
    NCO::SVector{NPIX,T}
    NC ::SVector{NPIX,T}
    xneq::SVector{N_neq,T}
end

const deff=1.0
const dust_to_gas_ratio = 1.0
const tdust = 15.
const tdust2 = tdust * 1e-2
const fa     = 1.0 / (1.0 + 1e4 * exp(-6e2 / tdust))

function get_kdust(temp)
    temp2  = temp * 1e-2

    stick  = 1.0 / (1.0 + 0.4 * (temp2 + tdust2)^0.5 + 0.2 * temp2 + 0.08 * temp2^2)
    ch7    = dust_to_gas_ratio * deff * 3e-17 * temp2^0.5
    Rdust  = ch7 * fa * stick
end

function calc_coeff!(coeff, par::Par, xelec)
    @unpack nH, temp, NH, NH2, NCO, NC, ξ, IUV, Zp = par

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
        coeff[iCion] *= fCselfshield(NC,NH2)
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

function init_abund(abund, Zp, xneq)
    #abund[dict["C+"]] = abC_s * Zp
    #abund[dict["O"]]  = abO_s * Zp
    calc_abund_derived(abund, Zp, xneq)
end

function calc_abund_derived(abund, Zp, xneq)
    if NONEQ
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

function solve_equilibrium_abundances(abund, dtime, par::Par)

    abund_dot = zeros(N_spec)
    abund_final = zeros(N_spec)
    abund_con = zeros(N_integrate)

    coeff = zeros(N_reac)
    reac_rates = zeros(N_reac)

    @unpack nH, Zp, xneq = par

    max_abund_inv = zeros(N_spec)
    get_max_abundance!(max_abund_inv, Zp)

    #we use closure s.t. large arrays like coeff & reac_rates can be updated in-place
    function calc_abund_dot_closure(abund_dot_int, abund_int, ppp, t)

        abund_dot .= 0.0 #in-place for better performance

        abund[idx_integrate] .= abund_int
        calc_abund_derived(abund, Zp, xneq)

        calc_coeff!(coeff, par, abund[ielec]) #zero allocation

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
    calc_abund_derived(abund_final, Zp, xneq)

    abund .= abund_final
    speciesoutofbound(abund, max_abund_inv, Inf) #fixing OFB error
    return sol.retcode, reac_rates
end


function speciesoutofbound(abund, max_abund_inv, ofbtol)
    for i in 1:N_spec
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

function get_max_abundance!(max_abund_inv, Zp)
    #update in-place
    for i in 1:N_spec
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

end #module
