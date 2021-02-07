#using DelimitedFiles
#using StaticArrays

ssCO = readdlm("CO_shielding_functions/shield.03.5.69-557-36.dat"); #CO linewidth = 0.3, T_ex = 5
const Nbin_nco = 46
const Nbin_nh2 = 41
const Nco = SVector{Nbin_nco,Float64}(ssCO[10:55])
const Nh2 = SVector{Nbin_nh2,Float64}(ssCO[58:98])
const istart = 100 #start line number for 12CO

function get_fco(is)
    a=[]
    for i in 0:4
        a = vcat(a, ssCO[is+i,:])
    end
    a = a[1:Nbin_nco]
end

ssCOtable = zeros(Nbin_nco, Nbin_nh2)
for j in 1:Nbin_nh2
    jj = istart+(j-1)*5
    ssCOtable[:, j] .= get_fco(jj)
end

#=
fCOss = interpolate((Nco, Nh2), ssCOtable, Gridded(Linear()));
function fCOselfshield(Nco,Nh2)
    Nh2 = Nh2 >= 1e23 ? 1e23 : Nh2
    Nco = Nco >= 1e19 ? 1e19 : Nco
    Nh2 = Nh2 <= 1e10 ? 1e10 : Nh2
    Nco = Nco <= 1e10 ? 1e10 : Nco
    return fCOss(Nco, Nh2)
end
=#

include("interpolation.jl")
const par_fCOss = InterpRange{Float64}(10., 19., 1/0.2, 15., 23., 1/0.2)
#const tblCO = SMatrix{Nbin_nco,Nbin_nh2}(ssCOtable) #size of SArray should not be larger than 100
const tblCO = Matrix{Float64}(ssCOtable)
function fCOselfshield(NCO,NH2)
    x, y = log10(NCO), log10(NH2)
    interp(x,y, tblCO, par_fCOss)
end


function fH2selfshield(Nh2)
    x = Nh2 * 2e-15
    b5 = 2.0
    return 0.965 / (1.0 + x/b5)^2 + 0.035 / sqrt(1 + x) * exp(-8.5e-4 * sqrt(1 + x))
end

function fCselfshield(NC,NH2)
    rH2 = 2.8e-22 * NH2
    return exp( -1.6e-17 * NC ) * exp(-rH2) / (1.0 + rH2)
end
