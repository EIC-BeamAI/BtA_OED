import SciBmad

println("Running: SciBmad public tracking integration"); flush(stdout)
@testset "SciBmad public tracking integration" begin
    beam = default_beam_moments()
    for apertures in (false, true)
        lattice = bta_lattice(; apertures)
        initial = reshape(beam.mean, 1, 6)
        reference = track_particles(initial; apertures)
        tracked = SciBmad.track(lattice; v0=initial, rf_on=false,
            use_explicit_SIMD=false, use_cpu_multithreading=false)
        @test tracked.bunch.coords.state == reference.coords.state
        @test tracked.bunch.coords.v ≈ reference.coords.v
        apertures || @test all(==(BeamTracking.STATE_ALIVE), reference.coords.state)
    end
end
