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

#read UMIST table
umist = readdlm("RATE12.dist.txt",':')

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
