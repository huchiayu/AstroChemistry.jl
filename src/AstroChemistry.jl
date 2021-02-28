module AstroChemistry

export solve_equilibrium_abundances, calc_abund_derived, init_abund, initialize_chemistry_network
export calc_coeff!
export Par, AbundTotal, Network
export fH2selfshield, fCOselfshield

using DelimitedFiles
using Printf
using StaticArrays
using DifferentialEquations
#using NLsolve
using Parameters
using Statistics

#using Sundials #for CVODE_BDF

include("OctTree.jl")
include("GadgetReader.jl")

include("shielding_functions.jl")
include("grain_recomb.jl")
include("constants.jl")
include("network.jl")
include("system_rate_eqs.jl")



end #module
