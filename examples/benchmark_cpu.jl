using BeamlineExperiments, BenchmarkTools

beam = default_beam_moments()
predict = prepared_prediction(prepare_bta())
cached = prepared_prediction(prepare_bta(); fixed_mean=beam.mean)
println("Fresh lattice:")
@btime predict_profiles($beam)
println("Prepared lattice:")
@btime $predict($beam)
println("Prepared transport tensors, beam covariance only:")
@btime $cached($beam)
