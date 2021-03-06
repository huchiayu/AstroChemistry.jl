#import Pkg; Pkg.add("HDF5"); Pkg.add("StaticArrays"); Pkg.add("PyPlot"); Pkg.add("DifferentialEquations"); Pkg.add("Parameters"); Pkg.add("Sundials");

using HDF5
using StaticArrays
using Statistics
using LinearAlgebra
using .Threads
using Serialization
using Random
using PyPlot
using Parameters
using BenchmarkTools

using AstroChemistry


const N = 3
const Nbin = 200
const xmin = 16.
const xmax = 23.


xl = "NH [cm^-2]"

const xmin = 16.
const xmax = 23.


function runPDR()

    T = Float64

    all_species::Vector{String} =
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

    net, dict = initialize_chemistry_network(all_species, grRec=false)
    @unpack iH2, iCO, iC, fac_H, fac_C, fac_O, charge = net
    abtot = AbundTotal(abC=1e-4, abO=3e-4, abSi=0.0)
    @unpack abC, abO, abSi = abtot


    nH = 1000.
    temp = 50.0
    ξ = 0.5e-16 #H2
    IUV = 10.0
    Zp = 1.0
    dtime = 1e10

    N_spec = length(all_species)
    @show N_spec
    ab_vs_x = Vector{Vector{T}}(undef, Nbin)
    for i in 1:Nbin
        ab_vs_x[i] = zeros(N_spec)
    end
    d_lnNH = (xmax-xmin) / (Nbin-1)
    xbin = 10 .^ (xmin .+ collect(0:Nbin-1) .* d_lnNH)
    NHbin = xbin
    dNH = NHbin[2:Nbin] .- NHbin[1:Nbin-1]

    NH2bin = zeros(Nbin)
    NCObin = zeros(Nbin)
    NCbin = zeros(Nbin)
    N_reac = length(net.alpha)
    reaction_rates=zeros(Nbin, N_reac)

    #@time for i in 1:1
    @time for i in eachindex(NHbin)
        xH2 = getindex.(ab_vs_x, iH2)
        xCO = getindex.(ab_vs_x, iCO)
        xC  = getindex.(ab_vs_x, iC )
        NH = NHbin[i]
        NH2 = NCO = NC = 1e10
        if i > 1
            NH2 = sum(xH2[1:i-1] .* dNH[1:i-1])
            NCO = sum(xCO[1:i-1] .* dNH[1:i-1])
            NC  = sum(xC[1:i-1]  .* dNH[1:i-1])
        end
        NH2bin[i] = NH2
        NCObin[i] = NCO
        NCbin[i] = NC
        #@show i, NH, NH2,NCO, NC
        #print(i," ")


        G = 2.8e-5 # for Zp = 1
        αG = IUV * AstroChemistry.kH2diss / AstroChemistry.kdust / nH * G
        #i==1 ? println("αG/2 = ", αG/2) : nothing

        xneq = SVector{0,T}([])
        #par = Par{1,0,T}(nH, temp, ξ, IUV, Zp, NH=SVector{1,T}(NH), NH2=SVector{1,T}(NH2), NCO=SVector{1,T}(NCO), NC=SVector{1,T}(NC), xneq=xneq)
        par = Par{1,0,T}(nH=nH, temp=temp, ξ=ξ, IUV=IUV, Zp=Zp,
        NH=[NH], NH2=[NH2], NCO=[NCO], NC=[NC], xneq=xneq)
        #NH=SVector{1,T}(NH), NH2=SVector{1,T}(NH2), NCO=SVector{1,T}(NCO), NC=SVector{1,T}(NC), xneq=xneq)

        retcode, reaction_rates[i,:] = solve_equilibrium_abundances(ab_vs_x[i], dtime, par, abtot, net)

        calc_abund_derived(ab_vs_x[i], xneq, abtot, net)
        sumH  = sum( ab_vs_x[i] .* fac_H )
        sumC  = sum( ab_vs_x[i] .* fac_C )
        sumO  = sum( ab_vs_x[i] .* fac_O )
        #sumSi = sum( ab_vs_x[i] .* fac_Si )
        #sumS  = sum( ab_vs_x[i] .* fac_S )
        sumelec = sum( ab_vs_x[i] .* charge )
        (1.0 ≈ sumH)  ? nothing : error("sumH = " , sumH)
        (abC  ≈ sumC)  ? nothing : error("sumC = " , sumC)
        (abO  ≈ sumO)  ? nothing : error("sumO = " , sumO)
        #(abSi ≈ sumSi) ? nothing : error("sumSi = ", sumSi)
        #(abS  ≈ sumS)  ? nothing : error("sumS = " , sumS)
        (1.0 ≈ 1.0 + sumelec) ? nothing : error("sumelec = ",sumelec)
    end

    return ab_vs_x, NHbin, NH2bin, NCObin, NCbin, reaction_rates, dict, net
end

ab_vs_x, N_H, N_H2, N_CO, N_C, rr, dict, net = runPDR();

xbin = 10 .^ (xmin .+ (xmax-xmin) .* collect(0:Nbin-1) ./ (Nbin-1));

clf()
fig, ax = PyPlot.subplots(1, 3, figsize=(15,5))

#abC_tot  = [sum(ab_vs_x[i] .* fac_C)  for i in 1:Nbin]
#abO_tot  = [sum(ab_vs_x[i] .* fac_O)  for i in 1:Nbin]
#abSi_tot = [sum(ab_vs_x[i] .* fac_Si) for i in 1:Nbin]
#abS_tot  = [sum(ab_vs_x[i] .* fac_S)  for i in 1:Nbin]

xminp, xmaxp = 20., 23.
ax[1].plot(xbin, 2 .*getindex.(ab_vs_x, dict["H2"]), "-"  , label="H2")
ax[1].plot(xbin, getindex.(ab_vs_x, dict["H"]), "--" , label="H")
ax[1].plot(xbin, getindex.(ab_vs_x, dict["e-"]), ":"  , label="e-")
ax[1].legend(loc="best", fontsize=9, ncol=1, frameon=false)
ax[1].set_xscale("log")
ax[1].set_yscale("log")
ax[1].axis([10.0^xminp, 10.0^xmaxp, 1e-6, 2])
ax[1].set_xlabel(xl)
ax[1].set_ylabel("x_i")
ax[1].grid(linestyle=":", linewidth=1)
ax2 = ax[1].twiny()
Av = N_H .* 5.35e-22
ax2.set_xlabel("Av")
ax2.set_xlim((10.0^xminp, 10.0^xmaxp).*5.35e-22)
ax2.set_xscale("log")
ax2.set_yscale("log")

xminp, xmaxp = 20., 23.
ax[2].plot(xbin, getindex.(ab_vs_x, dict["H+"]), "-" , label="H+")
ax[2].plot(xbin, getindex.(ab_vs_x, dict["H3+"]) , "--" , label="H3+")
ax[2].plot(xbin, getindex.(ab_vs_x, dict["He+"]) , "-." , label="He+")
ax[2].plot(xbin, getindex.(ab_vs_x, dict["O"])  , ":"  , label="O")
ax[2].plot(xbin, getindex.(ab_vs_x, dict["OH"]) , "-" , label="OH")
ax[2].plot(xbin, getindex.(ab_vs_x, dict["O2"]) , "--" , label="O2")
ax[2].plot(xbin, getindex.(ab_vs_x, dict["H2O"]), "-." , label="H2O")
ax[2].legend(loc="best", fontsize=9, ncol=2, frameon=false)
ax[2].set_xscale("log")
ax[2].set_yscale("log")
ax[2].axis([10.0^xminp, 10.0^xmaxp, 1e-11, 1e-3])
ax[2].set_xlabel(xl)
ax[2].grid(linestyle=":", linewidth=1)
ax2 = ax[2].twiny()
Av = N_H .* 5.35e-22
ax2.set_xlabel("Av")
ax2.set_xlim((10.0^xminp, 10.0^xmaxp).*5.35e-22)
ax2.set_xscale("log")
ax2.set_yscale("log")

ax[3].plot(xbin, getindex.(ab_vs_x, dict["C"])  , "-"  , label="C")
ax[3].plot(xbin, getindex.(ab_vs_x, dict["CO"]) , "--" , label="CO")
ax[3].plot(xbin, getindex.(ab_vs_x, dict["C+"]) , "-." , label="C+")
ax[3].plot(xbin, getindex.(ab_vs_x, dict["CH"]) , ":" , label="CH")
ax[3].plot(xbin, getindex.(ab_vs_x, dict["CH2"]) , ":" , label="CH2")
ax[3].plot(xbin, getindex.(ab_vs_x, dict["HCO+"]) , ":" , label="HCO+")
ax[3].legend(loc="upper left", fontsize=10, ncol=1, frameon=false)
ax[3].set_xscale("log")
ax[3].set_yscale("log")
ax[3].axis([10.0^xminp, 10.0^xmaxp, 1e-11, 1e-3])
ax[3].set_xlabel(xl)
ax[3].grid(linestyle=":", linewidth=1)
ax2 = ax[3].twiny()
Av = N_H .* 5.35e-22
ax2.set_xlabel("Av")
ax2.set_xlim((10.0^xminp, 10.0^xmaxp).*5.35e-22)
ax2.set_xscale("log")
ax2.set_yscale("log")

tight_layout()
savefig("PDRtest.png")
