"""
    UncoupledBeamState(; centroid, epsilon_x, beta_x, alpha_x, epsilon_y, beta_y, alpha_y)

Unknown transverse centroid `(x, px, y, py)` and uncoupled betatron second
moments represented by emittance/Twiss parameters. Positions are metres,
momenta are canonical normalized momenta, and geometric emittances are m⋅rad.
The example defaults match `default_beam_moments`, not a measured beam estimate.
Momentum spread, dispersion, mean delta and mean z are fixed simulation keywords.
"""
struct UncoupledBeamState{T<:Real}
    centroid::NTuple{4,T}
    epsilon_x::T
    beta_x::T
    alpha_x::T
    epsilon_y::T
    beta_y::T
    alpha_y::T
end

function UncoupledBeamState(; 
        centroid=(0.100092799, 0.0033, 0.0, 0.0),
        epsilon_x=1e-6, beta_x=12.6, alpha_x=1.8,
        epsilon_y=1e-6, beta_y=4.1, alpha_y=-0.65
    )
    length(centroid) == 4 || throw(DimensionMismatch("centroid must be (x, px, y, py)"))
    values = promote(centroid..., epsilon_x, beta_x, alpha_x, epsilon_y, beta_y, alpha_y, 0.0)
    all(isfinite, values) || throw(ArgumentError("beam state must be finite"))
    all(j -> values[j] > 0, (5, 6, 8, 9)) ||
        throw(ArgumentError("emittances and beta functions must be positive"))
    UncoupledBeamState(values[1:4], values[5:10]...)
end

function UncoupledBeamState(parameters::AbstractVector)
    length(parameters) == 10 || throw(DimensionMismatch("10 beam-state parameters required"))
    UncoupledBeamState(; centroid=parameters[1:4], epsilon_x=parameters[5],
        beta_x=parameters[6], alpha_x=parameters[7], epsilon_y=parameters[8],
        beta_y=parameters[9], alpha_y=parameters[10])
end

UncoupledBeamState(state::NamedTuple) = UncoupledBeamState(; state...)

"""Return `(x, px, y, py, epsilon_x, beta_x, alpha_x, epsilon_y, beta_y, alpha_y)`."""
uncoupled_beam_parameters(beam::UncoupledBeamState) =
    [beam.centroid..., beam.epsilon_x, beam.beta_x, beam.alpha_x,
     beam.epsilon_y, beam.beta_y, beam.alpha_y]

"""Column names for the 30-input multiwire Jacobian: 15 quad A, 5 steering A, 10 beam."""
bta_input_names() = vcat(
    [Symbol("quad_current_", name) for name in QUAD_NAMES],
    [Symbol("steering_current_", name) for name in CORRECTOR_NAMES],
    [:x, :px, :y, :py, :epsilon_x, :beta_x, :alpha_x, :epsilon_y, :beta_y, :alpha_y],
)

function pack_bta_inputs(quad_currents, steering_currents, beam::UncoupledBeamState)
    length(quad_currents) == 15 || throw(DimensionMismatch("15 quad currents required"))
    length(steering_currents) == 5 || throw(DimensionMismatch("5 steering currents required"))
    inputs = vcat(collect(quad_currents), collect(steering_currents), uncoupled_beam_parameters(beam))
    all(isfinite, inputs) || throw(ArgumentError("magnet currents and beam state must be finite"))
    inputs
end

pack_bta_inputs(quad_currents, steering_currents, beam::NamedTuple) =
    pack_bta_inputs(quad_currents, steering_currents, UncoupledBeamState(beam))

function _simulation_beam(
        beam::UncoupledBeamState;
        sigma_delta=1e-4, eta_x=2.7, eta_px=-0.42, delta0=0.0, z0=0.0
    )
    beam_moments_from_twiss(; mean=[beam.centroid..., z0, delta0],
        epsilon_x=beam.epsilon_x, beta_x=beam.beta_x, alpha_x=beam.alpha_x,
        epsilon_y=beam.epsilon_y, beta_y=beam.beta_y, alpha_y=beam.alpha_y,
        sigma_delta, eta_x, eta_px)
end

"""
    bta_measurements(inputs; ...)
    bta_measurements(quad_currents, steering_currents, beam_state; ...)

Deterministic loss-free multiwire summaries, `(mean_x, sigma_x, mean_y, sigma_y)`
in mm, one observation per row. `inputs` has the ordering from `bta_input_names()`.
This is the smooth function to compose with AD/OED. Its nested AD includes mixed
second tracking derivatives. Hard-aperture particle counting and noise are
provided separately by `simulate_bta`.
"""
function bta_measurements(
        inputs::AbstractVector;
        op=default_bta_operating_point(), observations=PROFILE_NAMES, moment_order=1,
        sigma_delta=1e-4, eta_x=2.7, eta_px=-0.42, delta0=0.0, z0=0.0,
        quad_factors=ones(15), quad_dx=zeros(15), quad_dy=zeros(15),
        calibration=DEFAULT_BTA_CALIBRATION
    )
    length(inputs) == 30 || throw(DimensionMismatch("30 simulation inputs required"))
    state = UncoupledBeamState(inputs[21:30])
    beam = _simulation_beam(state; sigma_delta, eta_x, eta_px, delta0, z0)
    setting = merge(op, (; quad_currents=inputs[1:15], corrector_currents=inputs[16:20]))
    predict_profiles(beam; op=setting, observations, moment_order, quad_factors,
                     quad_dx, quad_dy, calibration)
end

bta_measurements(quad_currents, steering_currents, beam_state; kwargs...) =
    bta_measurements(pack_bta_inputs(quad_currents, steering_currents, beam_state); kwargs...)

function _measurement_value_jacobian(inputs, options, ::Val{C}) where {C}
    f = p -> vec(bta_measurements(p; options...))
    config = ForwardDiff.JacobianConfig(f, inputs, ForwardDiff.Chunk{C}())
    ForwardDiff.jacobian(f, inputs, config)
end

function _summary_noise_scale(noise_std_mm, n)
    sigma = if noise_std_mm isa Real
        fill(Float64(noise_std_mm), n, 4)
    elseif noise_std_mm isa Union{Tuple,AbstractVector}
        length(noise_std_mm) == 4 || throw(DimensionMismatch("four channel noise scales required"))
        repeat(reshape(Float64.(collect(noise_std_mm)), 1, 4), n, 1)
    else
        size(noise_std_mm) == (n, 4) || throw(DimensionMismatch("noise must have shape ($n, 4)"))
        Float64.(noise_std_mm)
    end
    all(s -> isfinite(s) && s >= 0, sigma) ||
        throw(ArgumentError("noise standard deviations must be finite and nonnegative"))
    sigma
end

function _gaussian_particles(
        state::UncoupledBeamState, normals;
        sigma_delta, eta_x, eta_px, delta0, z0
    )
    # An explicit factor keeps the same base random draws across changes in the
    # beam state. It also supports zero momentum spread and a fixed z coordinate.
    L = zeros(6, 5)
    L[1, 1] = sqrt(state.epsilon_x * state.beta_x)
    L[2, 1] = -sqrt(state.epsilon_x / state.beta_x) * state.alpha_x
    L[2, 2] = sqrt(state.epsilon_x / state.beta_x)
    L[3, 3] = sqrt(state.epsilon_y * state.beta_y)
    L[4, 3] = -sqrt(state.epsilon_y / state.beta_y) * state.alpha_y
    L[4, 4] = sqrt(state.epsilon_y / state.beta_y)
    L[1, 5] = eta_x * sigma_delta
    L[2, 5] = eta_px * sigma_delta
    L[6, 5] = sigma_delta
    mean = [state.centroid..., z0, delta0]
    normals * transpose(L) .+ transpose(mean)
end

function _particle_summary(coordinates, alive)
    count_alive = count(alive)
    count_alive == 0 && return fill(NaN, 4)
    x = coordinates[alive, 1]
    y = coordinates[alive, 3]
    mean_x, mean_y = sum(x) / count_alive, sum(y) / count_alive
    # Population RMS of the equal-weight surviving particles, not sample stderr.
    [1000mean_x, 1000sqrt(sum(abs2, x .- mean_x) / count_alive),
     1000mean_y, 1000sqrt(sum(abs2, y .- mean_y) / count_alive)]
end

function _sampled_distribution(lattice, coordinates, observations)
    indices = _observation_indices(lattice, observations)
    bunch = Bunch(copy(coordinates); species=lattice.species_ref, p_over_q_ref=lattice.p_over_q_ref)
    n = size(coordinates, 1)
    summaries = fill(NaN, length(indices), 4)
    counts = zeros(Int, length(indices))
    first_loss = zeros(Int, n)
    for (j, ele) in enumerate(lattice.line)
        BeamTracking.track!(bunch, ele; rf_on=false, use_explicit_SIMD=false,
                            use_cpu_multithreading=false)
        alive = bunch.coords.state .== BeamTracking.STATE_ALIVE
        for k in 1:n
            !alive[k] && first_loss[k] == 0 && (first_loss[k] = j)
        end
        for row in findall(==(j), indices)
            counts[row] = count(alive)
            summaries[row, :] = _particle_summary(bunch.coords.v, alive)
        end
    end
    alive = BitVector(bunch.coords.state .== BeamTracking.STATE_ALIVE)
    final_distribution = (
        coordinates=copy(bunch.coords.v[alive, :]),
        particle_ids=findall(alive), alive,
        all_coordinates=copy(bunch.coords.v), state=copy(bunch.coords.state),
        first_loss_element=first_loss,
    )
    (; summaries, counts, final_distribution, lost=n - count(alive))
end

"""
    simulate_bta(quad_currents, steering_currents, beam_state=UncoupledBeamState(); ...)

One simulation call for design/control/inference. Inputs are 15 quadrupole
currents [A], five steering currents [A], and a transverse `UncoupledBeamState`
(or an equivalent keyword NamedTuple). The default tracks 1000 Gaussian particles.

Returns clean/noisy loss-free moment observations and their 16 × 30 Jacobian,
loss-aware particle observations at the harps, the surviving final distribution,
and fractional beam loss. Columns follow `bta_input_names()`, rows follow
`vec(observations.noiseless)` (monitor index varies fastest). Noise is independent
additive Gaussian summary noise in mm; the default 0.1 mm is illustrative.

`standard_particles`, if supplied, must be `nparticles × 5` standard-normal
variates, enabling common random numbers. Pass a seeded `rng` for reproducibility.
Noise is drawn before the particles so changing particle count does not change
the ideal noisy observation for a fixed RNG seed. `apertures=false` disables hard
apertures only for the particle path; the measurement/Jacobian path always ignores
loss. Set `compute_jacobian=false` to skip the Jacobian when only simulation is needed.

Lost particles' `final_distribution.all_coordinates` are last tracked states,
not end-of-line coordinates. Only `final_distribution.coordinates` contains
survivors at the end. If no particle reaches a harp, its particle summary is NaN.
"""
function simulate_bta(quad_currents, steering_currents,
        beam_state::Union{UncoupledBeamState,NamedTuple}=UncoupledBeamState();
        op=default_bta_operating_point(), observations=PROFILE_NAMES,
        noise_std_mm=0.1, rng::AbstractRNG=Random.default_rng(), nparticles::Integer=1000,
        standard_particles=nothing, compute_jacobian::Bool=true, chunksize::Integer=4,
        moment_order=1, sigma_delta=1e-4, eta_x=2.7, eta_px=-0.42, delta0=0.0, z0=0.0,
        quad_factors=ones(15), quad_dx=zeros(15), quad_dy=zeros(15),
        calibration=DEFAULT_BTA_CALIBRATION, apertures::Bool=true
    )
    nparticles > 0 || throw(ArgumentError("nparticles must be positive"))
    chunksize > 0 || throw(ArgumentError("chunksize must be positive"))
    state = beam_state isa NamedTuple ? UncoupledBeamState(beam_state) : beam_state
    inputs = Float64.(pack_bta_inputs(quad_currents, steering_currents, state))
    sigma = _summary_noise_scale(noise_std_mm, length(observations))
    all(isfinite, (sigma_delta, eta_x, eta_px, delta0, z0)) && sigma_delta >= 0 ||
        throw(ArgumentError("fixed beam settings must be finite; sigma_delta must be nonnegative"))
    if !isnothing(standard_particles)
        size(standard_particles) == (nparticles, 5) ||
            throw(DimensionMismatch("standard_particles must be ($nparticles, 5)"))
        all(isfinite, standard_particles) || throw(ArgumentError("standard particles must be finite"))
    end
    options = (; op, observations, moment_order, sigma_delta, eta_x, eta_px,
                delta0, z0, quad_factors, quad_dx, quad_dy, calibration)
    clean = bta_measurements(inputs; options...)
    jacobian = compute_jacobian ?
        _measurement_value_jacobian(inputs, options, Val(min(chunksize, 30))) : nothing
    perturbation = sigma .* randn(rng, size(clean)...)
    noisy = clean + perturbation
    normals = isnothing(standard_particles) ? randn(rng, nparticles, 5) : Float64.(standard_particles)
    # The particle path is deliberately primal; differentiate bta_measurements.
    primal_state = UncoupledBeamState(inputs[21:30])
    initial_particles = _gaussian_particles(primal_state, normals;
        sigma_delta, eta_x, eta_px, delta0, z0)
    setting = merge(op, (; quad_currents=inputs[1:15], corrector_currents=inputs[16:20]))
    lattice = bta_lattice(setting; apertures, quad_factors, quad_dx, quad_dy, calibration)
    tracked = _sampled_distribution(lattice, initial_particles, observations)
    (
        observations=(names=Tuple(Symbol.(observations)),
            columns=(:mean_x, :sigma_x, :mean_y, :sigma_y), units=:mm,
            noiseless=clean, noisy, noise_std_mm=sigma,
            particle_noiseless=tracked.summaries,
            particle_noisy=tracked.summaries + perturbation,
            particle_counts=tracked.counts),
        jacobian, jacobian_input_names=bta_input_names(),
        jacobian_reference=:lossless_moment_observations, inputs,
        final_distribution=tracked.final_distribution,
        fractional_beam_loss=tracked.lost / nparticles,
        n_particles=nparticles, n_lost=tracked.lost,
        initial_distribution=initial_particles,
    )
end
