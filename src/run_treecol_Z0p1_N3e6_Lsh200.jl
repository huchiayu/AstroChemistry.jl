#species_noneq = []
#species_noneq = ["H2"]
species_noneq = ["H2", "H+"]
const N_neq = length(species_noneq)
const NONEQ = N_neq > 0

const grRec = true

const file_path = "/ptmp/huchiayu/snapshots/tallbox/SFSNPI_N3e6_gS10H250dS40H250_soft0p3_from200Myr_SFcut_Z0p1_Lsh200"
const Zp = 0.1
const ShieldingLength = 0.2
include("treecol.jl")

#const Zp = parse(T, ARGS[2])

@show file_path
@show Zp

#snaps = collect(740:-10:150)
snaps = collect(550:10:1000)
for i in snaps
    abund_all, ga, X, NH, NH2, NCO, tree = solve_chem_all_particles(i, file_path, Zp);
end
0



