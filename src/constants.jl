const XHe = 0.1
const omega = 0.5;
const SEC_PER_YEAR = 31556952.0

#const abC_s = 2.9e-4
#const abO_s = 4.9e-4
#const abN_s = 6.8e-4
#const abSi_s = 3.2e-5
#const abC_s = 1.4e-4
#const abO_s = 3.2e-4
#const abS_s = 8.0e-6
#const abN_s = 7.5e-5
#const abMg_s = 3.5e-5
#const abFe_s = 2.8e-5
#const abF_s = 2.e-8
#const abNa_s = 2.e-8
const abC_s = 1.4e-4
const abO_s = 3.2e-4
#const abSi_s = 0.0
#const abSi_s = 1.5e-5
const abSi_s = 1.7e-6 #depleted abundance
const abS_s = 0.0
#const abS_s = 1e-6
#const abN_s = 7.6e-5
#const abFe_s = 1.6e-6
#const abMg_s = 1.4e-5
#const abNa_s = 1.3e-6
const abN_s = 0.0
const abFe_s = 0.0
const abMg_s = 0.0
const abNa_s = 0.0
const abF_s = 0.0



#always present
const iH = dict["H"]
const ielec = dict["e-"]
const iH2 = dict["H2"]
const iHp = dict["H+"]

#optional elements
const iHe = any(all_species.=="He") ? dict["He"] : 0
const iC  = any(all_species.=="C" ) ? dict["C" ] : 0
const iO  = any(all_species.=="O" ) ? dict["O" ] : 0
const iSi = any(all_species.=="Si") ? dict["Si"] : 0
const iS  = any(all_species.=="S" ) ? dict["S" ] : 0
const iN  = any(all_species.=="N" ) ? dict["N" ] : 0
const iMg = any(all_species.=="Mg" ) ? dict["Mg" ] : 0
const iFe = any(all_species.=="Fe" ) ? dict["Fe" ] : 0

const iCO  = any(all_species.=="CO" ) ? dict["CO" ] : 0
const iHep = any(all_species.=="He+") ? dict["He+"] : 0
const iCp = any(all_species.=="C+") ? dict["C+"] : 0
const iSip = any(all_species.=="Si+") ? dict["Si+"] : 0


const INTEGRATE = true  # 1 = integration, 0 = root-finding

const ξ0 = 1.36e-17
const kdust = 3e-17
#const kH2diss = 5.68e-11 * 0.52
const kH2diss = 5.68e-11
const kCOdiss = 2.43e-10 * 0.48
const gamma_H2 = 3.85
const gamma_CO = 3.51

const iCion = 0
const SupTh = false
const grRec = false
#const grRec = true
const ssC = false

#=
species_noneq = []
#species_noneq = ["H2"]
#species_noneq = ["H2", "H+"]

const N_neq = length(species_noneq)
const NONEQ = N_neq > 0
=#
