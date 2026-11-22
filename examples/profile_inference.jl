using BeamlineExperiments, ForwardDiff

op = default_bta_operating_point()
beam = default_beam_moments()
data = reshape(predict_profiles(beam; op), 1, 4, 4) .+ 0.05
# Differentiate a synthetic harp likelihood at a fixed lattice.
decode = beam_decoder((:x, :px, :δ0, :η_x, :σ_δ))
experiments = prepare_bta_batch([op])
loglikelihood = bta_loglikelihood(experiments, data, decode; noise_std_mm=0.1)
θ = [beam.mean[1:2]..., 0.0, 2.7, 1e-4]
@show loglikelihood(θ) ForwardDiff.gradient(loglikelihood, θ)
