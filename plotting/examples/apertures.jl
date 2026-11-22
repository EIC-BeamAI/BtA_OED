using GLMakie
include(joinpath(@__DIR__, "..", "BeamlineExperimentsPlotting.jl"))
using .BeamlineExperimentsPlotting
GLMakie.activate!()
input = isempty(ARGS) ? default_lattice() : ARGS[1]
view = interactive_apertures(input)
# `julia -i` keeps the REPL available; ordinary script execution waits for closing.
if !isinteractive()
    wait(display(view.figure))
end
