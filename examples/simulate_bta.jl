using BtaUQ
using Random

op = default_bta_operating_point()
beam = UncoupledBeamState(
    centroid=(0.100092799, 0.0033, 0.0, 0.0),
    epsilon_x=1e-6, beta_x=12.6, alpha_x=1.8,
    epsilon_y=1e-6, beta_y=4.1, alpha_y=-0.65,
)

result = simulate_bta(op.quad_currents, op.corrector_currents, beam;
    rng=MersenneTwister(42), nparticles=1000,
    # Illustrative independent summary uncertainties, all in mm:
    noise_std_mm=(0.1, 0.05, 0.1, 0.05))

println("Multiwires: ", result.observations.names)
println("Columns: ", result.observations.columns, " [mm]")
println("Clean, loss-free moment observations:")
display(result.observations.noiseless)
println("Noisy observations:")
display(result.observations.noisy)
println("Jacobian shape: ", size(result.jacobian))
println("Jacobian columns: ", result.jacobian_input_names)
println("Particles reaching each harp: ", result.observations.particle_counts)
println("Loss-aware particle summaries:")
display(result.observations.particle_noiseless)
println("Fractional beam loss: ", result.fractional_beam_loss,
        " (", result.n_lost, "/", result.n_particles, ")")
println("Surviving end-of-line distribution: ", size(result.final_distribution.coordinates))

# Jacobian is d vec(clean observations) / d [quad currents; steer currents; beam state].
# Each column's row order is MW006/MW060/MW125/MW166 for x, then sigma_x, y, sigma_y.
quad_jacobian = result.jacobian[:, 1:15]
steering_jacobian = result.jacobian[:, 16:20]
beam_jacobian = result.jacobian[:, 21:30]
