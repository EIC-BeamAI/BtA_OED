using CairoMakie
include(joinpath(@__DIR__, "..", "BeamlineExperimentsPlotting.jl"))
using .BeamlineExperimentsPlotting
CairoMakie.activate!()
input = isempty(ARGS) ? default_lattice() : ARGS[1]
output = length(ARGS) >= 2 ? ARGS[2] : joinpath(@__DIR__,"..","output","bta-apertures.png")
mkpath(dirname(abspath(output)))
view = interactive_apertures(input; show=false)
save(output,view.figure; px_per_unit=1.5)
println("Saved ", abspath(output))
