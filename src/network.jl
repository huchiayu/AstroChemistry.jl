@with_kw struct Network{N_spec, N_reac, N_integrate, N_neq, T}

    species_noneq::Vector{String}
    grRec::Bool
    SupTh::Bool
    ssC  ::Bool

    iCion     :: Int32
    iH2formGr :: Int32
    iH2diss   :: Int32
    iCOdiss   :: Int32
    iHpRecGr  :: Int32
    iHepRecGr :: Int32
    iCpRecGr  :: Int32
    iSipRecGr :: Int32

    fac_H ::SVector{N_spec,Int32}
    fac_He::SVector{N_spec,Int32}
    fac_C ::SVector{N_spec,Int32}
    fac_O ::SVector{N_spec,Int32}
    fac_Si::SVector{N_spec,Int32}
    fac_S ::SVector{N_spec,Int32}
    fac_N ::SVector{N_spec,Int32}
    fac_Mg::SVector{N_spec,Int32}
    fac_Fe::SVector{N_spec,Int32}
    charge::SVector{N_spec,Int32}

    typ  :: SVector{N_reac,Int32}
    ir1  :: SVector{N_reac,Int32}
    ir2  :: SVector{N_reac,Int32}
    ip1  :: SVector{N_reac,Int32}
    ip2  :: SVector{N_reac,Int32}
    ip3  :: SVector{N_reac,Int32}
    ip4  :: SVector{N_reac,Int32}
    pro3 :: Bool
    pro4 :: Bool

    alpha :: SVector{N_reac,T}
    beta  :: SVector{N_reac,T}
    gamma :: SVector{N_reac,T}

    #always present
    iH    :: Int32
    ielec :: Int32
    iH2   :: Int32
    iHp   :: Int32

    #optional elements
    iHe   :: Int32
    iC    :: Int32
    iO    :: Int32
    iSi   :: Int32
    iS    :: Int32
    iN    :: Int32
    iMg   :: Int32
    iFe   :: Int32
    iCO   :: Int32
    iHep  :: Int32
    iCp   :: Int32
    iSip  :: Int32

    ineq          :: SVector{N_neq,Int32}
    idx_integrate :: SVector{N_integrate,Int32}
end


function read_umist_line(umist_line)
    num = umist_line[1]
    typ = umist_line[2]
    R::Vector{String} = [umist_line[3], umist_line[4]]
    P::Vector{String} = [umist_line[5], umist_line[6], umist_line[7], umist_line[8]]
    #alpha, beta, gamma, Tl, Tu, ST, ACC =
    #    umist_line[10], umist_line[11], umist_line[12],
    #    umist_line[13], umist_line[14], umist_line[15], umist_line[16]
    alpha, beta, gamma =
        umist_line[10], umist_line[11], umist_line[12]
    return num, typ, R, P, alpha, beta, gamma
end

function print_reaction(num, typ, R, P, alpha, beta, gamma)
    @printf "(%4s)" num
    #@printf "[%3s]" typ

    @printf " [α,β,γ = "
    @printf "%10.2e," alpha
    @printf "%10.2e," beta
    @printf "%10.2e]" gamma

    @printf "%8s" R[1]
    @printf "     + %8s" R[2]
    if typ == "GR"
        @printf "(gr)-->%8s" P[1]
    else
        @printf "    -->%8s" P[1]
    end
    if(P[2]!="")
        @printf "     + %8s" P[2]
    end
    if(P[3]!="")
        @printf "     + %8s" P[3]
    end
    if(P[4]!="")
        @printf "     + %8s" P[4]
    end
    print("\n")
end

function print_reaction(num, typ, R, P)
    @printf "%7s" R[1]
    @printf "     + %7s" R[2]
    if typ == "GR"
        @printf " (gr)-->%7s" P[1]
    else
        @printf "     -->%7s" P[1]
    end
    if(P[2]!="")
        @printf "     + %7s" P[2]
    end
    if(P[3]!="")
        @printf "     + %7s" P[3]
    end
    if(P[4]!="")
        @printf "     + %7s" P[4]
    end
    print("\n")
end

function search_species_reactants(network, s1::String)
    reac=[]
    for i in 1:size(network, 1)
        n, t, R, P, a, b, c = read_umist_line(network[i,:])
        if any(s1 .== R)
            print_reaction(n, t, R, P, a, b, c)
            push!(reac, n)
        end
    end
    reac
end

function search_species_reactants(network, s1::String, s2::String)
    reac=[]
    for i in 1:size(network, 1)
        n, t, R, P, a, b, c = read_umist_line(network[i,:])
        if any(s1 .== R) && any(s2 .== R)
            @printf "(i=%.4d) " i
            print_reaction(n, t, R, P, a, b, c)
            push!(reac, n)
        end
    end
    reac
end

function search_species_products(network, s1::String)
    reac=[]
    for i in 1:size(network, 1)
        n, t, R, P, a, b, c = read_umist_line(network[i,:])
        if any(s1 .== P)
            print_reaction(n, t, R, P, a, b, c)
            push!(reac, n)
        end
    end
    reac
end

function search_species_products(network, s1::String, s2::String)
    reac=[]
    for i in 1:size(network, 1)
        n, t, R, P, a, b, c = read_umist_line(network[i,:])
        if any(s1 .== P) && any(s2 .== P)
            print_reaction(n, t, R, P, a, b, c)
            push!(reac, n)
        end
    end
    reac
end

function print_all_reactions(network)
    for i in 1:size(network, 1)
        n, t, R, P, a, b, c = read_umist_line(network[i,:])
        #@printf "(%3d)  " i
        print_reaction(n, t, R, P, a, b, c)
    end
end

function get_all_species_in_network(network)
    species::Vector{String} = []
    for i in 1:size(network, 1)
        n, t, R, P, a, b, c = read_umist_line(network[i,:])
        #@printf "(%3d)  " i
        #print_reaction(n, t, R, P, a, b, c)
        [push!(species, R[j]) for j in eachindex(R)]
        [push!(species, P[k]) for k in eachindex(P)]
    end
    species = unique(species)
    filter!(s -> (s != "")&&(s != "PHOTON")&&(s != "CRPHOT")&&(s != "CRP"), species)
end

function breakdown_species(s::String)
    charge::Int32 = 0
    elements = String[]
    counts   = Int32[]
    mult::Int32 = 1
    skipflag::Bool = false

    if s == "e-"
        return -1, elements, counts
    end

    for i in length(s):-1:1
        if skipflag
            #@show skipflag
            skipflag = false
            continue
        end
        char = s[i]
        #@show char
        isuppercase(char)
        if char == '+'
            charge += 1
        elseif char == '-'
            charge -= 1
        elseif isletter(char)
            elem = ""
            if islowercase(char)
                elem = s[i-1:i]
                skipflag = true
            else
                elem = string(char)
            end
            push!(elements, elem)
            push!(counts, mult)
            mult = 1
        else
            mult = parse(Int32, char)
        end
    end
    return charge, elements, counts
end

#s="HSi3Mg2+"
#breakdown_species(s)


function get_unique_elements(elements::Vector{Vector{String}})
    unique_elements = String[]
    for i in 1:length(elements)
        for j in eachindex(elements[i])
            push!(unique_elements, elements[i][j])
        end
    end
    return unique!(unique_elements)
end

function get_idx_integrate(elements::Vector{Vector{String}}, species_noneq::Vector{String}, dict::Dict)
    N_spec = length(elements)
    unique_elements = get_unique_elements(elements)
    #const N_uni_elem = length(unique_elements_a)
    #const uni_elem = SVector{N_uni_elem}(unique_elements_a)
    #N_conserve = length(unique_elements) + 1  #+1 for electron

    #species_derived = ["e-", "H", "He", "C", "O", "Si", "S"]
    species_derived = unique_elements
    push!(species_derived, "e-")
    #@show species_derived

    N_neq = length(species_noneq)
    if N_neq > 0
        for i in 1:N_neq
            push!(species_derived, species_noneq[i])
        end
    end
    ineq = SVector{N_neq}([dict[species_noneq[i]] for i in 1:N_neq])
    #@show species_derived
    #@show ineq
    #@show all_species
    idx_derived = sort([dict[species_derived[i]] for i in eachindex(species_derived)])

    idx_integrate_a = filter!(s -> all(s .!= idx_derived), collect(1:N_spec))
    return idx_integrate_a, species_derived
end
function construct_index_arrays(reduced_network::Matrix{Any}, dict::Dict)

    N_reac = size(reduced_network, 1)
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

    iCion     = 0
    iH2formGr = 0
    iH2diss   = 0
    iCOdiss   = 0
    iHpRecGr  = 0
    iHepRecGr = 0
    iCpRecGr  = 0
    iSipRecGr = 0

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
        if t=="GR" && any(R .== "H") && any(P .== "H2") #H2 photodissociation
            iH2formGr = i
        end
        if any(R .== "H2") && any(R .== "PHOTON") && any(P .== "H") #H2 photodissociation
            iH2diss = i
        end
        if any(R .== "CO") && any(R .== "PHOTON") && any(P .== "C") && any(P .== "O") #CO photodissociation
            iCOdiss = i
        end
        if t=="GR" && any(R .== "H+") && any(R .== "e-") && any(P .== "H") && any(P .== "PHOTON") #H+ recomb. on grains
            iHpRecGr = i
        end
        if t=="GR" && any(R .== "He+") && any(R .== "e-") && any(P .== "He") && any(P .== "PHOTON") #He+ recomb. on grains
            iHepRecGr = i
        end
        if t=="GR" && any(R .== "C+") && any(R .== "e-") && any(P .== "C") && any(P .== "PHOTON") #C+ recomb. on grains
            iCpRecGr = i
        end
        if t=="GR" && any(R .== "Si+") && any(R .== "e-") && any(P .== "Si") && any(P .== "PHOTON") #H+ recomb. on grains
            iSipRecGr = i
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
        #if SupTh == true && t == "IN"
        #    if any(R .== "C+") && any(R .== "H2") && any(P .== "CH+") && any(P .== "H")
        #        typ_a[i] = 5
        #        @show i,R,P
        #    end
        #end

        ir1_a[i], ir2_a[i] = dict[R[1]], dict[R[2]]
        ip1_a[i], ip2_a[i], ip3_a[i], ip4_a[i] = dict[P[1]], dict[P[2]], dict[P[3]], dict[P[4]]
    end
    return typ_a, ir1_a, ir2_a, ip1_a, ip2_a, ip3_a, ip4_a, alpha_a, beta_a, gamma_a,
    iCion, iH2formGr, iH2diss, iCOdiss, iHpRecGr, iHepRecGr, iCpRecGr, iSipRecGr
end
function setup_conservation_arrays(elem::String, elements::Vector{Vector{String}}, counts::Vector{Vector{Int32}})
    N_spec = length(elements)
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

function initialize_chemistry_network(all_species::Vector{String};
    species_noneq::Vector{String} = String[],
    grRec::Bool = true,
    SupTh::Bool = false,
    ssC  ::Bool = false,
    T::DataType = Float64
    )
    N_spec = length(all_species)

    all_species_2 = vcat(all_species , ["", "PHOTON", "CRP", "CRPHOT"]);

    #read UMIST table
    umist = readdlm((@__DIR__)*"/RATE12.dist.txt",':')

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
    more_reactions = readdlm((@__DIR__)*"/user-rates.txt",':')
    reac_num=[]
    for i in 1:size(more_reactions,1)
        #println(i)
        n, t, R, P, a, b, c = read_umist_line(more_reactions[i,:])
        if any(R[1].==all_species_2) && any(R[2].==all_species_2) && any(P[1].==all_species_2) && any(P[2].==all_species_2) && any(P[3].==all_species_2) && any(P[4].==all_species_2)
            #print_reaction(n, t, R, P, a, b, c)
            push!(reac_num,n)
        end
    end
    #@show reac_num
    more_reactions = more_reactions[reac_num,:]
    reduced_network = vcat(reduced_network, more_reactions);
    print_all_reactions(reduced_network)

    N_reac = size(reduced_network, 1)
    #@show N_reac, N_spec

    dict = Dict(all_species .=> collect(1:N_spec) );
    dict[""] = dict["PHOTON"] = dict["CRP"] = dict["CRPHOT"] = 0

    charge = @. Int32(getindex(breakdown_species(all_species) , 1))
    elements = @. getindex(breakdown_species(all_species) , 2)
    counts = @. getindex(breakdown_species(all_species) , 3)

    idx_integrate, species_derived = get_idx_integrate(elements, species_noneq, dict)
    N_integrate = length(idx_integrate)

    N_neq = length(species_noneq)
    ineq = SVector{N_neq,Int32}([dict[species_noneq[i]] for i in 1:N_neq])

    typ, ir1, ir2, ip1, ip2, ip3, ip4, alpha, beta, gamma,
    iCion, iH2formGr, iH2diss, iCOdiss, iHpRecGr, iHepRecGr, iCpRecGr, iSipRecGr = construct_index_arrays(reduced_network, dict)

    net = Network{N_spec, N_reac, N_integrate, N_neq, Float64}(
    grRec = grRec,
    SupTh = SupTh,
    ssC   = ssC,
    iCion     = iCion    ,
    iH2formGr = iH2formGr,
    iH2diss   = iH2diss  ,
    iCOdiss   = iCOdiss  ,
    iHpRecGr  = iHpRecGr ,
    iHepRecGr = iHepRecGr,
    iCpRecGr  = iCpRecGr ,
    iSipRecGr = iSipRecGr,
    species_noneq = species_noneq,
    fac_H  = setup_conservation_arrays("H" , elements, counts),
    fac_He = setup_conservation_arrays("He", elements, counts),
    fac_C  = setup_conservation_arrays("C" , elements, counts),
    fac_O  = setup_conservation_arrays("O" , elements, counts),
    fac_Si = setup_conservation_arrays("Si", elements, counts),
    fac_S  = setup_conservation_arrays("S" , elements, counts),
    fac_N  = setup_conservation_arrays("N" , elements, counts),
    fac_Mg = setup_conservation_arrays("Mg", elements, counts),
    fac_Fe = setup_conservation_arrays("Fe", elements, counts),
    charge = charge,
    typ    = typ  ,
    ir1    = ir1  ,
    ir2    = ir2  ,
    ip1    = ip1  ,
    ip2    = ip2  ,
    ip3    = ip3  ,
    ip4    = ip4  ,
    pro3   = any(ip3 .!= 0),
    pro4   = any(ip4 .!= 0),
    alpha  = alpha,
    beta   = beta ,
    gamma  = gamma,
    #always present
    iH = dict["H"],
    ielec = dict["e-"],
    iH2 = dict["H2"],
    iHp = dict["H+"],
    iHe = any(all_species.=="He") ? dict["He"] : 0,
    iC  = any(all_species.=="C" ) ? dict["C" ] : 0,
    iO  = any(all_species.=="O" ) ? dict["O" ] : 0,
    iSi = any(all_species.=="Si") ? dict["Si"] : 0,
    iS  = any(all_species.=="S" ) ? dict["S" ] : 0,
    iN  = any(all_species.=="N" ) ? dict["N" ] : 0,
    iMg = any(all_species.=="Mg" ) ? dict["Mg" ] : 0,
    iFe = any(all_species.=="Fe" ) ? dict["Fe" ] : 0,
    iCO  = any(all_species.=="CO" ) ? dict["CO" ] : 0,
    iHep = any(all_species.=="He+") ? dict["He+"] : 0,
    iCp = any(all_species.=="C+") ? dict["C+"] : 0,
    iSip = any(all_species.=="Si+") ? dict["Si+"] : 0,
    #N_neq = N_neq,
    #N_integrate = N_integrate,
    idx_integrate = idx_integrate,
    ineq = ineq
    )
    return net, dict
end


@with_kw struct AbundTotal{T}
    abC :: T = 1.4e-4
    abO :: T = 3.2e-4
    # abSi = 0.0
    # abSi = 1.5e-5
    abSi :: T = 1.7e-6 #depleted abundance
    # abS = 1e-6
    # abN = 7.6e-5
    # abFe = 1.6e-6
    # abMg = 1.4e-5
    abS :: T = 0.0
    abN :: T = 0.0
    abFe :: T = 0.0
    abMg :: T = 0.0
end

import Base.+
+(a::AbundTotal{T}, b::AbundTotal{T}) where {T<:Real} =
AbundTotal{T}(
a.abC  + b.abC ,
a.abO  + b.abO ,
a.abSi + b.abSi,
a.abS  + b.abS ,
a.abN  + b.abN ,
a.abFe + b.abFe,
a.abMg + b.abMg)

import Base.-
-(a::AbundTotal{T}, b::AbundTotal{T}) where {T<:Real} =
AbundTotal{T}(
a.abC  - b.abC ,
a.abO  - b.abO ,
a.abSi - b.abSi,
a.abS  - b.abS ,
a.abN  - b.abN ,
a.abFe - b.abFe,
a.abMg - b.abMg)

import Base.*
*(a::AbundTotal, c::T) where {T<:Real} = AbundTotal(c*a.abC, c*a.abO, c*a.abSi, c*a.abS, c*a.abN, c*a.abFe, c*a.abMg)
*(c::T, a::AbundTotal) where {T<:Real} = *(a::AbundTotal, c::T)
import Base./
/(a::AbundTotal, c::T) where {T<:Real} = AbundTotal(a.abC/c, a.abO/c, a.abSi/c, a.abS/c, a.abN/c, a.abFe/c, a.abMg/c)


#=
function test()
    all_species =
    String[
    "H", "H-", "H2", "H+", "H2+", "H3+", "e-", "He", "He+", "HeH+",
    "C", "C+", "CO", "HCO+",
    "O", "O+", "OH", "OH+", "H2O+", "H3O+", "H2O",
    "O2", "CO+", "O2+",
    "CH2", "CH2+", "CH", "CH+", "CH3+",
    "Si+", "Si"
    #"SiO", "SiO+", "SiOH+",
    #"S", "S+",
    ]

    species_noneq=String[]
    #net::Network{31,286,25,0,Float64}, dict = initialize_chemistry_network(all_species, species_noneq)
    net, dict = initialize_chemistry_network(all_species, species_noneq)

    #@code_warntype initialize_chemistry_network(all_species, species_noneq)
    @unpack iH2, iCO, iC, fac_H, fac_C, fac_O, charge = net
    abtot = AbundTotal()
    @unpack abC_s, abO_s, abSi_s = abtot

    net
end
@code_warntype test()
#n = test()
=#

###### specific reactions
#=
reactions = Int32[]

rea_H = [6155, 75, 4921, 6093, 489, 731, 730, 2604, 1278, 1277, 2895, 733, 406, 408, 434, 6162, 492]
rea_He = [734, 3431, 491, 458, 2677, 2610, 2678, 3065, 6094,6166, 3517, 6173, 967]
rea_O = [2686, 2951, 2689, 2662, 1296, 1295, 1293, 1294, 1266, 1264, 1265, 1427, 5643, 5199,
1646, 4029, 408, 1613, 2789, 2767, 2904, 378, 5196, 1643, 1642, 5705, 948, 405, 1421, 5412]
rea_CH = [2852, 2647, 6095, 2648, 1175, 1176, 1177, 1158, 1161, 1160, 1159, 3055, 5370, 5374]
rea_CO = [5199, 1646, 488, 2654, 1349, 5310, 5231, 5232, 64, 1349, 823, 6158, 707, 194, 876]
#rea_photo = [5997, 5915, 5989, 5887, 5888, 5902]
rea_suprathermal = [2618, 2189]

reactions = vcat(rea_H, rea_He, rea_O, rea_CH, rea_CO, rea_suprathermal)
unique!(reactions);

reduced_network = umist[reactions,:]
const all_species = get_all_species_in_network(reduced_network);

#add photoprocesses automatically
charge_a = @. Int32(getindex(breakdown_species(all_species) , 1));
photo_reactions = []
for i in 1:size(all_species,1)
    #if charge_a[i] <= 0
    if true
        r = search_species_reactants(umist, "PHOTON", all_species[i])
        if all_species[i] != "CO" #will do CO shielding manually
            global photo_reactions = vcat(r, photo_reactions)
        end
    end
end

#add photoprocesses to the network
all_reactions = vcat(reactions, photo_reactions)
unique!(all_reactions);
reduced_network = umist[all_reactions,:];

const N_reac = size(reduced_network, 1)
const all_species = get_all_species_in_network(reduced_network);
const N_spec = size(all_species,1)
#println(all_species)
#print_all_reactions(reduced_network)
#@show N_reac, N_spec
=#

######
