using Test, CairoMakie, PALSParserJ
include(joinpath(@__DIR__, "..", "plotting", "BeamlineExperimentsPlotting.jl"))
using .BeamlineExperimentsPlotting
CairoMakie.activate!()
@testset "Artifact-backed BeamlineExperiments plotting" begin
    @test !haskey(ENV,"PALS_PARSER_CPP_LIB")
    @test !haskey(ENV,"PALS_PARSER_CPP_DIR")
    @test evaluate_pals_expression("2+3*4") == 14
    @test occursin("/artifacts/",replace(PALSParserJ.libparser(),'\\'=>'/'))
    view = interactive_apertures(; show=false)
    @test length(view.records) == 80
    @test length(view.tab.name) == 109
    @test maximum(view.tab.s) ≈ 64.70348665922872
    @test view.slider.value[] == 30
    original = copy(view.data[].highedge)
    actual_profile = copy(view.data[].profile_hi)
    set_close_to!(view.slider,60)
    @test !isequal(original,view.data[].highedge)
    @test isequal(actual_profile,view.data[].profile_hi)
    # The opening doubles when the exaggeration is doubled.
    oldgeom = BeamlineExperimentsPlotting.ribbon_geometry(view.tab,view.records,30)
    newgeom = view.data[]
    width(d) = hypot((d.highedge[1] - d.lowedge[1])...)
    @test width(newgeom) ≈ 2width(oldgeom) rtol=1e-5
    for scale in (1,100)
        set_close_to!(view.slider,scale)
        @test view.slider.value[] == scale
    end
    view.reset.clicks[] += 1
    @test_throws ErrorException interactive_apertures(;scale=0,show=false)
    mktempdir() do dir
        file = joinpath(dir,"apertures.png")
        save(file,view.figure)
        @test filesize(file)>10_000
        fp = floor_view()
        save(joinpath(dir,"floor.png"),fp.figure)
        @test filesize(joinpath(dir,"floor.png"))>1_000
    end
end
