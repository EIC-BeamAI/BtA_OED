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
    @test length(view.tab.name) == 110
    @test maximum(view.tab.s) ≈ 64.70348665922872
    @test view.slider.value[] == 30
    # The terminal bore uses the translated exit frame rather than the
    # extrapolated entrance axis.
    true_scale = BeamlineExperimentsPlotting.ribbon_geometry(view.tab, view.records, 1)
    i = findlast(==("patch_l20_3"), view.tab.name)
    theta = view.tab.theta[i]
    origin = [view.tab.z[i], view.tab.x[i]] +
        0.20875625 .* [cos(theta), sin(theta)] +
        0.08362315 .* [-sin(theta), cos(theta)]
    theta_out = theta - 9.7789e-5
    midpoint = collect((true_scale.patch_openings[end-1] + true_scale.patch_openings[end]) / 2)
    @test midpoint ≈ origin - 0.01853565 .* [-sin(theta_out), cos(theta_out)] atol=1e-5 rtol=0
    @test collect(last(true_scale.centerline)) ≈ origin atol=1e-5 rtol=0
    # Project the two outward edges onto the same transverse axis.
    outward_edge(d) = collect(d.patch_openings[end-1])
    preceding_edge = collect(true_scale.patch_openings[end-3])
    transverse = [-sin(theta), cos(theta)]
    step = sum((outward_edge(true_scale) - preceding_edge) .* transverse)
    @test step ≈ -0.005022 atol=1e-5 rtol=0
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
    exit_midpoint(d) = collect((d.patch_openings[end-1] + d.patch_openings[end]) / 2)
    exit_width(d) = hypot((d.patch_openings[end] - d.patch_openings[end-1])...)
    @test exit_midpoint(newgeom) ≈ midpoint atol=1e-5 rtol=0
    @test exit_width(newgeom) ≈ 60exit_width(true_scale) rtol=1e-5
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
