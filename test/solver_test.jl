push!(LOAD_PATH, pwd())
using StaticArrays
using Statistics
using LinearAlgebra
using .Threads
using Serialization
using Random

#species_noneq = []
#species_noneq = ["H2"]
species_noneq = ["H2", "H+"]
const N_neq = length(species_noneq)
const NONEQ = N_neq > 0

#include("AstroChemistry.jl")
using AstroChemistry

const T = Float64
const N = 3


dt = 1e9

abund = zeros(N_spec)
nH, temp, ξ, IUV, Zp = 4219., 1127., 1.15e-16, 0.888, 1.0
NH_eff, NH2_eff, NCO_eff = 1.157e22, 5.08e21, -0.0
#@show nH, temp, ξ, IUV, Zp, NH_eff, NH2_eff, NCO_eff
npix = 1
#xneq = SVector{N_neq,T}(0.4973, 4.1089e-13)
xneq = SVector{N_neq,T}(0.4997, 4.1089e-13)
#xneq = SVector{N_neq,T}(0.4997, 3.711e-9)
#xneq = SVector{N_neq,T}(0.4997)
#xneq = SVector{N_neq,T}()

abund[dict["CO"]] = abC_s * Zp
#abund[dict["C+"]] = abC_s * Zp - abund[dict["CO"]]
#abund[dict["H2"]] = xneq[1]
#abund[dict["H+"]] = xneq[2]
init_abund(abund, Zp, xneq)

par = Par{npix,T}(nH, temp, ξ, IUV, Zp,
    SVector{npix,T}(NH_eff),
    SVector{npix,T}(NH2_eff),
    SVector{npix,T}(NCO_eff),
    SVector{npix,T}(0), xneq)

@time retcode, rr = solve_equilibrium_abundances(abund, dt, par)
