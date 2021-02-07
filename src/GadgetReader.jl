#module GadgetReader
using HDF5
export read_snap

function read_snap(filename)

    T=Float64
    header::Dict = h5readattr(filename, "/Header")
    boxsize::T = header["BoxSize"]
    time::T    = header["Time"]

    N_gas::Int64 = header["NumPart_ThisFile"][1]

    pos_gas::Matrix{T} = h5read(filename, "PartType0/Coordinates");
    vel_gas::Matrix{T} = h5read(filename, "PartType0/Velocities");
    rho::Vector{T}     = h5read(filename, "PartType0/Density");
    u::Vector{T}       = h5read(filename, "PartType0/InternalEnergy");
    m_gas::Vector{T}   = h5read(filename, "PartType0/Masses");
    hsml::Vector{T}    = h5read(filename, "PartType0/SmoothingLength");
    #scal::Vector{T}       = h5read(filename, "PartType0/PassiveScalarField");

    id_gas::Vector{Int64} = h5read(filename, "PartType0/ParticleIDs");
    abund::Matrix{T} = h5read(filename, "PartType0/ChemicalAbundancesSG");
    fH2ss::Vector{T}    = h5read(filename, "PartType0/ShieldingFactorH2");
    fH2dust::Vector{T}    = h5read(filename, "PartType0/ShieldingFactorDust");

    col_tot::Matrix{T}    = h5read(filename, "PartType0/TreecolColumnDensitiesAll");
    col_H2::Matrix{T}    = h5read(filename, "PartType0/TreecolColumnDensitiesH2");
    col_CO::Matrix{T}    = h5read(filename, "PartType0/TreecolColumnDensitiesCO");

    Tdust::Vector{T}    = h5read(filename, "PartType0/DustTemperature");

    N_star::Int64 = header["NumPart_ThisFile"][5]

    pos_star::Matrix{T} = h5read(filename, "PartType4/Coordinates");
    vel_star::Matrix{T} = h5read(filename, "PartType4/Velocities");
    m_star::Vector{T}   = h5read(filename, "PartType4/Masses");
    sftime::Vector{T}      = h5read(filename, "PartType4/StellarFormationTime");
    id_star::Vector{Int64} = h5read(filename, "PartType4/ParticleIDs");

    return N_gas, pos_gas, vel_gas, rho, u, m_gas, hsml, id_gas,
        abund, fH2ss, fH2dust, col_tot, col_H2, col_CO, Tdust,
        N_star, pos_star, vel_star, m_star, sftime, id_star,
        boxsize, time
end

#end
