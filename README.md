# AstroChemistry

```AstroChemistry.jl``` is a package that solves a user-defined chemistry network for astrophysical applications (in particular for the interstellar medium). It can be used to post-process snapshots of hydrodynamical simulations with a built-in Octree that calculates the column densities for shielding against the FUV radiation. The chemistry coefficients are based on the [UMIST](http://udfa.ajmarkwick.net/index.php) database and the CO shielding is based on [Visser et al. 2009](https://home.strw.leidenuniv.nl/~ewine/photo/CO_photodissociation.html). The code is multithreading parallel using Julia's built-in functionality. 

# Dependencies
```AstroChemistry.jl``` makes use of the following Julia packages:
- [StaticArrays.jl](https://github.com/JuliaArrays/StaticArrays.jl)
- [Parameters.jl](https://github.com/mauro3/Parameters.jl)
- [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl)
- [Healpix.jl](https://github.com/ziotom78/Healpix.jl)

# Installation
Open the REPL (type ```julia``` in the terminal) and install the relevant packages:
```
import Pkg
Pkg.add("StaticArrays")
Pkg.add("Parameters")
Pkg.add("DifferentialEquations")
Pkg.add("Healpix")
Pkg.add("AstroChemistry")
```
# Quick Start
```
using AstroChemistry

#specify the chemical species in the network
all_species = ["H", "H-", "H2", "H+", "H2+", "e-"]

#initialize the network
net, dict = initialize_chemistry_network(all_species, grRec=false)

#specify the physical parameters
nH, temp, ξ, IUV, Zp = 100., 20., 1e-16, 0., 1.
par = Par{1,0,Float64}(nH=nH, temp=temp, ξ=ξ, IUV=IUV, Zp=Zp)

#initialize the total (gas-phase) elemental abundance
abtot = AbundTotal() * Zp

#initialize the chemical abundance array
N_spec = length(all_species)
abund = zeros(N_spec)
init_abund(abund, [], abtot, net)

#solve the network
dt = 1e10 #yr
retcode, rr = solve_equilibrium_abundances(abund, dt, par, abtot, net)

#print the final abundances of H+, H and H2
@show retcode, abund[dict["H+"]], abund[dict["H"]], abund[dict["H2"]]
```

# Examples

- 1D PDR calculation: see [examples/PDRtest.jl](https://github.com/huchiayu/AstroChemistry.jl/blob/main/examples/PDRtest.jl) for reproducing the "F1" test problem described in the PDR code comparison project [Röllig et al. 2017](https://www.aanda.org/component/article?access=bibcode&bibcode=&bibcode=2007A%2526A...467..187RFUL).

- post-processing: see [examples/chem_octree.jl](https://github.com/huchiayu/AstroChemistry.jl/blob/main/examples/chem_octree.jl) for how to read a [Gadget](https://wwwmpa.mpa-garching.mpg.de/gadget4/) snapshot (format 3) and calculate the chemical abundances in post-processing.

# Gallery
- The output image of [examples/PDRtest.jl](https://github.com/huchiayu/AstroChemistry.jl/blob/main/examples/PDRtest.jl).
![PDRtest](https://user-images.githubusercontent.com/23061774/109493462-d1467d80-7a8c-11eb-94e4-3f03252bbf2c.png)

- Post-processing results from a hydrodynamical simulation of the interstellar medium in the solar-neighborhood ([Hu et al. 2021](https://arxiv.org/abs/2103.03889)).
![load_mol_maps_neq_620](https://user-images.githubusercontent.com/23061774/109887710-82facf80-7c82-11eb-8753-085ae225e497.png)

# Author
Chia-Yu Hu @ Max Planck Institute for Extraterrestrial Physics 
(cyhu.astro@gmail.com)
