push!(LOAD_PATH, "../src")

using HDF5
using StaticArrays
using Statistics
using LinearAlgebra
using .Threads
using Serialization
using Random
using Healpix
using Parameters

#include("GadgetReader.jl")
#include("OctTree.jl")
using GadgetReader
using OctTree

using AstroChemistry

const N = 3


function vec2svec(vec::Matrix{T}) where {T}
    svec = [SVector{3,T}(vec[:,i]) for i in 1:size(vec,2)]
end
function mat2smat(mat::Array{T,3}) where {T}
    smat = [SMatrix{3,3,T}(mat[:,:,i]) for i in 1:size(mat,3)]
end



const ANGLE = 0.7
const ShieldingLength = 1e10

const XH = 0.71
const BOLTZMANN=1.3806e-16
const PROTONMASS=1.6726e-24
const GRAVCON=6.67e-8
const UnitMass_in_g = 1.989e43
const UnitLength_in_cm    = 3.085678e21
const UnitTime_in_s = 3.08568e+16

const UnitDensity_in_cgs = UnitMass_in_g / UnitLength_in_cm^3
const UnitDensity_in_pccm = UnitDensity_in_cgs/PROTONMASS
const Year_in_s = 31556926.

const fac_col = (UnitMass_in_g/UnitLength_in_cm^2)/PROTONMASS

const facNHtoAv = 5.35e-22
const gamma_eff = 3.51

function solve_chem_all_particles(file_path, Zp)

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
    N_spec = length(all_species)
    N_neq = 0

    net, dict = initialize_chemistry_network(all_species)
    @unpack iH2, iCO, iC, fac_H, fac_C, fac_O, charge = net
    abtot = AbundTotal()
    @unpack abC_s, abO_s, abSi_s = abtot


    fname = file_path * "/snapshot.hdf5"
    Npart, pos, vel, rho, u, mass, hsml, id_gas,
    boxsize, time = read_snap(fname);
    println("boxsize = ", boxsize)

    boxsizes = SVector(boxsize, boxsize, boxsize)
    center = 0.5 .* boxsizes
    topnode_length = boxsizes

    X = vec2svec(pos);


    ξ = 1.3e-16  #H2
    IUV = 1.0

    mu = 2.3 #mean molecular weight
    nH = rho .* (XH * UnitDensity_in_cgs / PROTONMASS)
    temp = u ./ (1.5 * BOLTZMANN /mu /PROTONMASS / 1e10)

    mass_H2 = zeros(Npart)
    mass_CO = zeros(Npart)

    #NH = NH2 = NCO = NC = 0.0 #don't declear here!! otherwise they won't be thread-safe

    abund_all = [zeros(N_spec) for _ in 1:Npart]


    for i in 1:Npart
        xneq = SVector{N_neq,T}()
        #make sure the abC is consistent between simulations & postprocessing
        abund_all[i][dict["C+"]] = abC_s * Zp
        #initial guess for H2 & H+ (useful for calcultating steady state H2)
        #abund_all[i][dict["H2"]] = 0.49
        init_abund(abund_all[i], Zp, xneq, abtot, net)
    end

    NH_eff = zeros(Npart)
    NH2_eff = zeros(Npart)
    NCO_eff = zeros(Npart)
    ga_out = Vector{TreeGather{Float64}}(undef,Npart)
    tree_out = nothing
    status = zeros(Int, Npart)


    Nstep = 3
    time_max = 1e9 #unit = yr
    dt = time_max / Nstep
    NCpix = zeros(NPIX)

    for j in 1:Nstep
        @. mass_H2 = mass * getindex(abund_all,iH2) * XH * 2.0
        @. mass_CO = mass * getindex(abund_all,iCO) * XH * 28.0

        println("iteration ", j)

        println("buildtree...")
        @time tree = buildtree(X, hsml, mass, mass_H2, mass_CO, center, topnode_length);

        println("loop over particles and solve the chemistry network...")
        ichem = findall((nH.>0.1).&(temp.<5e3))
        println(length(ichem), " cold and dense particles found...")
        @time @threads for i in ichem #better load balance
            if i%1000 == 0
                print("i=", i, " ")
            end
            ga = TreeGather{T}()
            treewalk(ga,X[i],tree,ANGLE,ShieldingLength,boxsizes)
            ga_out[i] = ga
            NH = ga.column_all .* (fac_col*XH)
            NH2 = ga.column_H2 .* (fac_col/2)
            NCO = ga.column_CO .* (fac_col/28)
            NC = 0.0
            facNHtoAv_eff = gamma_eff*facNHtoAv*Zp
            NH_eff[i]  = -log(mean(exp.(-NH .*facNHtoAv_eff))) / facNHtoAv_eff
            NH2_eff[i] = -log(mean(exp.(-NH2.*facNHtoAv_eff))) / facNHtoAv_eff
            NCO_eff[i] = -log(mean(exp.(-NCO.*facNHtoAv_eff))) / facNHtoAv_eff

            xneq = SVector{N_neq,T}()
            par = Par{NPIX,N_neq,T}(nH[i], temp[i], ξ, IUV, Zp,
                SVector{NPIX,T}(NH),
                SVector{NPIX,T}(NH2),
                SVector{NPIX,T}(NCO),
                SVector{NPIX,T}(NCpix), xneq)
            #@show NH,NH2,NCO,NCpix,xneq
            #dt = 1e9 / (Zp * nH[i]) #in years
            retcode, _  = solve_equilibrium_abundances(abund_all[i], dt, par, abtot, net)
            if retcode != :Success
                status[i] = 1  #1 means failed
                @show retcode
                @show i, nH[i], temp[i], ξ, IUV, Zp
                @show NH_eff[i], NH2_eff[i], NCO_eff[i], xneq
                #@show NH, NH2, NCO
            end
        end
        println("chemsitry done!")
        @show sum(status)

        xCO=getindex.(abund_all, dict["CO"])
        println("CO mass = ",  1e10*sum(mass.*xCO) * XH * 28 )
        xH2=getindex.(abund_all, dict["H2"])
        println("H2 mass = ",  1e10*sum(mass.*xH2) * XH * 2  )

        #this may seem redundant as it'll be computed in the next iteration
        #however, treewalk is so fast (compared to chemistry) that it doesn't hurt
        #and we get the updated column densities for little overhead
        println("get column densities for the outputs...")
        @. mass_H2 = mass * getindex(abund_all,iH2) * XH * 2.0
        @. mass_CO = mass * getindex(abund_all,iCO) * XH * 28.0
        println("buildtree...")
        @time tree = buildtree(X, hsml, mass, mass_H2, mass_CO, center, topnode_length);
        println("loop over particles and get column densities...")
        @time @threads for i in eachindex(NH_eff)
            ga = TreeGather{T}()
            treewalk(ga,X[i],tree,ANGLE,ShieldingLength,boxsizes)
            ga_out[i] = ga
            NH = ga.column_all .* (fac_col*XH)
            NH2 = ga.column_H2 .* (fac_col/2)
            NCO = ga.column_CO .* (fac_col/28)
            facNHtoAv_eff = gamma_eff*facNHtoAv*Zp
            NH_eff[i]  = -log(mean(exp.(-NH .*facNHtoAv_eff))) / facNHtoAv_eff
            NH2_eff[i] = -log(mean(exp.(-NH2.*facNHtoAv_eff))) / facNHtoAv_eff
            NCO_eff[i] = -log(mean(exp.(-NCO.*facNHtoAv_eff))) / facNHtoAv_eff
        end
        println("done!")
        tree_out = tree


        # write to file
        abund_all_arr = zeros(N_spec,Npart)
        for i in eachindex(abund_all)
            abund_all_arr[:,i] = abund_all[i]
        end
        #T = Float64
        fnamebase = "/data-chem-"
        if net.grRec == true
            fnamebase *= "GrRec-"
        end
        fname = file_path * fnamebase * "-" * string(j) *".hdf5"
        fid=h5open(fname,"w")
        grp_head = create_group(fid,"Header");
        attributes(fid["Header"])["all_species"] = all_species
        attributes(fid["Header"])["time"] = time
        grp_part = create_group(fid,"Chemistry");
        h5write(fname, "Gas/Position"        , pos)
        h5write(fname, "Gas/Density"         , nH)
        h5write(fname, "Gas/Temperature"     , temp)
        h5write(fname, "Gas/ID"              , id_gas)
        h5write(fname, "Chemistry/Abundances"     , abund_all_arr)
        h5write(fname, "Chemistry/NH_eff"         , NH_eff)
        h5write(fname, "Chemistry/NH2_eff"        , NH2_eff)
        h5write(fname, "Chemistry/NCO_eff"        , NCO_eff)
        h5write(fname, "Chemistry/Status"         , status)
        close(fid)
    end #iteration loop
    return pos, nH, temp, id_gas, abund_all, ga_out, X, NH_eff, NH2_eff, NCO_eff, tree_out
end

file_path = "./"
Zp = 1.0
pos, nH, temp, id_gas, abund_all, ga_out, X, NH_eff, NH2_eff, NCO_eff, tree_out = solve_chem_all_particles(file_path, Zp)
