push!(LOAD_PATH, pwd())
using HDF5
using StaticArrays
using Statistics
using LinearAlgebra
using .Threads
using Serialization
using Random

#include("OctTree.jl")
#using OctTree
#using GadgetReader
#include("AstroChemistry.jl")
using AstroChemistry

using Healpix


const T = Float64
const N = 3


function vec2svec(vec::Matrix{T}) where {T}
    svec = [SVector{3,T}(vec[:,i]) for i in 1:size(vec,2)]
end
function mat2smat(mat::Array{T,3}) where {T}
    smat = [SMatrix{3,3,T}(mat[:,:,i]) for i in 1:size(mat,3)]
end



const ANGLE = 0.7
#const ShieldingLength = 0.5

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

function solve_chem_all_particles(i, file_path, Zp)
    println("=============== snapshot ", i, " ===============")

    snap = ""
    if i < 10
        snap = "00" * string(i)
    elseif i < 100
        snap = "0" * string(i)
    else
        snap = string(i)
    end

    #X = [SVector{N}( rand(MersenneTwister(i),N) ) for i in 1:Npart-1]
    #push!(X,SVector(0.,0.,0.))
    #push!(X,SVector(0.5,0.5,0.5))

    fname = file_path * "/snap_" * snap * ".hdf5"
    Npart, pos, vel, rho, u, mass, hsml, id_gas,
        abund, fH2, fdust, col_tot, col_H2, col_CO, Tdust,
        N_star, pos_star, vel_star, m_star, sftime, id_star,
        boxsize, time = read_snap(fname);
    println("boxsize = ", boxsize)

    boxsizes = SVector(boxsize, boxsize, 20.0)
    center = SVector(0.5*boxsize, 0.5*boxsize, 0.0)
    topnode_length = SVector(10., 10., 10.)


    #down-sampling
    #=
    Ns::Int = 100
    rng = MersenneTwister(1114);
    ridx = randperm(rng, Npart);
    Npart = div(Npart,Ns)
    ridx = ridx[1:Npart]
    pos=pos[:,ridx]
    vel=vel[:,ridx]
    rho=rho[ridx]
    u=u[ridx]
    mass=mass[ridx]
    hsml=hsml[ridx]
    id_gas=id_gas[ridx]
    abund=abund[:,ridx]
    mass .*= Ns
    =#

    X = vec2svec(pos);
    #X .*= BOXSIZE
    #push!(X,SVector(0.,0.,0.) .+ 0.01)

    dt::T = 3e-2   #30Myr
    idx = @. ( time - sftime < dt )
    SFR = 1e10 * sum(abs.(m_star[idx])) / (dt*1e9)
    facSFR = SFR / 2.5e-3
    #facSFR = 1.0

    ξ = 1.3e-16 * facSFR #H2
    IUV = 1.0 * (facSFR > 2e-3 ? facSFR : 2e-3) #cosmic UV background
    @show facSFR, SFR, ξ, IUV

    mu = 2.3 #mean molecular weight
    nH = rho .* (XH * UnitDensity_in_cgs / PROTONMASS)
    temp = u ./ (1.5 * BOLTZMANN /mu /PROTONMASS / 1e10)

    mass_H2 = zeros(Npart)
    mass_CO = zeros(Npart)

    #NH = NH2 = NCO = NC = 0.0 #don't declear here!! otherwise they won't be thread-safe

    abund_all = [zeros(N_spec) for _ in 1:Npart]

    if NONEQ==false
        println("steady state: set a ceiling xH2 for the initial guess")
        abund[1, abund[1,:] .> 0.45] .= 0.45  #improve stability
    end

    for i in 1:Npart
        #abund = abund_all[:,1]
        #xneq = SVector{1,T}([abund[1,i]])
        if NONEQ
            xneq = SVector{N_neq,T}([abund[1,i], abund[2,i]])
        else
            xneq = SVector{N_neq,T}()
        end
        #use C+ & CO from simulations as the initial guess
        #abund_all[i][dict["CO"]] = abund[3,i] * (abC_s/3.01e-4)
        #abund_all[i][dict["C+"]] = abC_s * Zp - abund_all[i][dict["CO"]]
        #make sure the abC is consistent between simulations & postprocessing
        abund_all[i][dict["C+"]] = abC_s * Zp
        #initial guess for H2 & H+ (useful for calcultating steady state H2)
        abund_all[i][dict["H2"]] = abund[1,i]
        abund_all[i][dict["H+"]] = abund[2,i]
        init_abund(abund_all[i], Zp, xneq)
    end

    #hsml = ones(Npart) .* (0.5 * boxsizes[1])
    NH_eff = zeros(Npart)
    NH2_eff = zeros(Npart)
    NCO_eff = zeros(Npart)
    #ga_out = TreeGather{T}()
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
        #@time tree = buildtree(X, hsml, mass, mass_H2, mass_CO, boxsizes);
        @time tree = buildtree(X, hsml, mass, mass_H2, mass_CO, center, topnode_length);

        println("loop over particles and solve the chemistry network...")
        #@time for i in 273432:273432
        ichem = findall((nH.>1).&(temp.<3e3))
        println(length(ichem), " cold and dense particles found...")
        @time @threads for i in ichem #better load balance
        #@time for i in ichem #better load balance
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

            if NONEQ
                xneq = SVector{N_neq,T}([abund[1,i], abund[2,i]])
            else
                xneq = SVector{N_neq,T}()
            end
            temp[i] = temp[i] > 3.0 ? temp[i] : 3.0
            #nH[i] = nH[i] < 3e5 ? nH[i] : 3e5
            par = Par{NPIX,T}(nH[i], temp[i], ξ, IUV, Zp,
                SVector{NPIX,T}(NH),
                SVector{NPIX,T}(NH2),
                SVector{NPIX,T}(NCO),
                SVector{NPIX,T}(NCpix), xneq)
            #@show NH,NH2,NCO,NCpix,xneq
            #dt = 1e9 / (Zp * nH[i]) #in years
            retcode, _  = solve_equilibrium_abundances(abund_all[i], dt, par)
            if retcode != :Success
                status[i] = 1  #1 means failed
                @show retcode
                @show i, nH[i], temp[i], ξ, IUV, Zp
                @show NH_eff[i], NH2_eff[i], NCO_eff[i], xneq
                #@show NH, NH2, NCO
                #@show abund[1,i],abund[2,i]
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
        #@time tree = buildtree(X, hsml, mass, mass_H2, mass_CO, boxsizes);
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
        if NONEQ
            fnamebase = "/chem-neqH2Hp-noCOic-TF3-"
        else
            fnamebase = "/chem-eqH2Hp-noCOic-TF3-"
        end
        if grRec
            fnamebase *= "GrRec-"
        end
        fname = file_path * fnamebase * snap * "-" * string(j) *".hdf5"
        fid=h5open(fname,"w")
        grp_head = g_create(fid,"Header");
        attrs(fid["Header"])["all_species"] = all_species
        attrs(fid["Header"])["time"] = time
        attrs(fid["Header"])["facSFR"] = facSFR
        grp_part = g_create(fid,"Chemistry");
        h5write(fname, "Chemistry/Abundances"     , abund_all_arr)
        h5write(fname, "Chemistry/ID"             , id_gas)
        h5write(fname, "Chemistry/NH_eff"         , NH_eff)
        h5write(fname, "Chemistry/NH2_eff"        , NH2_eff)
        h5write(fname, "Chemistry/NCO_eff"        , NCO_eff)
        h5write(fname, "Chemistry/Status"         , status)
        close(fid)
    end #iteration loop
    return abund_all, ga_out, X, NH_eff, NH2_eff, NCO_eff, tree_out
end

#=
ix,iy = 1,2
clf()
fig = figure("",figsize=(8,8))
#plot_quadtree(tree, ix, iy)
scatter(getindex.(X,ix), getindex.(X,iy), c="blue", marker="o", s=0.3)
plot_treewalk(ga,ix,iy)
p = X[Npart]
plot(p[ix],p[iy],"o")
gcf()
#savefig("particles.png")
=#


#===
p = boxsizes .* 0.5

m = Map{Float64, RingOrder}(NSIDE);
m.pixels[:] = log10.(column_all)
#Plots.plot(m, clim=(3e4,8e4))
Plots.plot(m)

column_true = zeros(12*NSIDE^2)
@time for i in 1:length(X)
    dx = OctTree.nearest.(X[i] - p, boxsizes)
    ipix = OctTree.vec2pix(dx)
    #ga.column_all[ipix] += node.mass
    area = 4*pi/(NSIDE^2*12.) * sum(dx.^2)
    column_true[ipix] += mass[i] / area
end

m.pixels[:] = log10.(column_true)
Plots.plot(m)
===#
