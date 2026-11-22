# Requires a desktop OpenGL context; separate from the headless suite.
using Test, GLMakie
include(joinpath(@__DIR__, "..", "plotting", "BeamlineExperimentsPlotting.jl"))
using .BeamlineExperimentsPlotting
GLMakie.activate!()
@testset "Interactive aperture OpenGL rendering" begin
    view = interactive_apertures(; show=false)
    screen = GLMakie.Screen(; visible=false)
    try
        display(screen,view.figure)
        before = copy(GLMakie.Makie.colorbuffer(screen))
        @test length(unique(before))>1
        set_close_to!(view.slider,60)
        after = GLMakie.Makie.colorbuffer(screen)
        @test count(before .!= after)>0
        @test view.slider.value[]==60
    finally
        close(screen)
    end
end
