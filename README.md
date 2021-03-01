# AstroChemistry

```AstroChemistry.jl``` is a package that solves a user-difined chemistry network in the astrophysical context, in particular for the interstellar medium. It can be used to post-process snapshots of hydrodynamical simulations with an built-in Octree that calculates the column densities for shielding against the FUV radiation.

# Examples

- 1D PDR calculation: see [examples/PDRtest.jl](https://github.com/huchiayu/AstroChemistry.jl/blob/main/examples/PDRtest.jl) for reproducing the "F1" test problem described in the PDR code comparison project [Röllig et al. 2017](https://www.aanda.org/component/article?access=bibcode&bibcode=&bibcode=2007A%2526A...467..187RFUL).
- post-processing: see [examples/chem_octree.jl](https://github.com/huchiayu/AstroChemistry.jl/blob/main/examples/chem_octree.jl) for how to read a [Gadget](https://wwwmpa.mpa-garching.mpg.de/gadget4/) snapshot (format 3) and calculate the chemical abundances in post-processing.

