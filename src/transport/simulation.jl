"""Entrance centroid `(x, px, y, py, delta)` and uncoupled betatron optics.

Horizontal dispersion and RMS relative momentum spread determine the
momentum-dependent part of the covariance. Longitudinal position is fixed
separately in transport.
"""
struct UncoupledBeamState{T<:Real}
    centroid::NTuple{5,T}
    epsilon_x::T
    beta_x::T
    alpha_x::T
    epsilon_y::T
    beta_y::T
    alpha_y::T
    eta_x::T
    eta_px::T
    sigma_delta::T
end

const BEAM_PARAMETER_NAMES = (:x, :px, :y, :py, :δ0, :β_x, :α_x, :β_y, :α_y,
    :η_x, :η_px, :ε_x, :ε_y, :σ_δ)

function UncoupledBeamState(; centroid=(0.1013, 0.00174, 0.0, 0.0, 0.0),
                            epsilon_x=1e-6, beta_x=12.6, alpha_x=1.8,
                            epsilon_y=1e-6, beta_y=4.1, alpha_y=-0.65,
                            eta_x=2.7, eta_px=-0.42, sigma_delta=1e-4)
    length(centroid) == 5 || throw(DimensionMismatch("centroid must be (x, px, y, py, delta)"))
    values = promote(centroid..., epsilon_x, beta_x, alpha_x, epsilon_y, beta_y,
                     alpha_y, eta_x, eta_px, sigma_delta)
    all(isfinite, values) || throw(ArgumentError("beam state must be finite"))
    all(j -> values[j] > 0, (6, 7, 9, 10)) ||
        throw(ArgumentError("emittances and beta functions must be positive"))
    sigma_delta >= 0 || throw(ArgumentError("momentum spread must be nonnegative"))
    UncoupledBeamState(values[1:5], values[6:end]...)
end

function UncoupledBeamState(parameters::AbstractVector)
    length(parameters) == length(BEAM_PARAMETER_NAMES) ||
        throw(DimensionMismatch("expected $(length(BEAM_PARAMETER_NAMES)) beam parameters"))
    UncoupledBeamState(; centroid=parameters[1:5],
        beta_x=parameters[6], alpha_x=parameters[7], beta_y=parameters[8],
        alpha_y=parameters[9], eta_x=parameters[10], eta_px=parameters[11],
        epsilon_x=parameters[12], epsilon_y=parameters[13], sigma_delta=parameters[14])
end

UncoupledBeamState(state::NamedTuple) = UncoupledBeamState(; state...)

"""Pack the entrance state in BEAM_PARAMETER_NAMES order."""
uncoupled_beam_parameters(beam::UncoupledBeamState) =
    [beam.centroid..., beam.beta_x, beam.alpha_x, beam.beta_y, beam.alpha_y,
     beam.eta_x, beam.eta_px, beam.epsilon_x, beam.epsilon_y, beam.sigma_delta]

"""Column names for the multiwire Jacobian: quad currents, corrector currents, and beam."""
bta_input_names() = vcat(
    [Symbol("quad_current_", name) for name in QUAD_NAMES],
    [Symbol("corrector_current_", name) for name in CORRECTOR_NAMES],
    collect(BEAM_PARAMETER_NAMES),
)

function pack_bta_inputs(quad_currents, corrector_currents, beam::UncoupledBeamState)
    length(quad_currents) == 15 || throw(DimensionMismatch("15 quad currents required"))
    length(corrector_currents) == 5 || throw(DimensionMismatch("5 corrector currents required"))
    inputs = vcat(collect(quad_currents), collect(corrector_currents), uncoupled_beam_parameters(beam))
    all(isfinite, inputs) || throw(ArgumentError("magnet currents and beam state must be finite"))
    inputs
end

pack_bta_inputs(quad_currents, corrector_currents, beam::NamedTuple) =
    pack_bta_inputs(quad_currents, corrector_currents, UncoupledBeamState(beam))

_beam_optics(beam::UncoupledBeamState) =
    (; beam.epsilon_x, beam.beta_x, beam.alpha_x, beam.epsilon_y, beam.beta_y,
       beam.alpha_y, beam.eta_x, beam.eta_px, beam.sigma_delta)

_canonical_mean(beam::UncoupledBeamState, z0) =
    [beam.centroid[1:4]..., z0, beam.centroid[5]]

beam_moments_from_twiss(beam::UncoupledBeamState; z0=0.0, mean=_canonical_mean(beam, z0)) =
    beam_moments_from_twiss(; mean, _beam_optics(beam)...)

"""
    bta_measurements(inputs; ...)
    bta_measurements(quad_currents, corrector_currents, beam_state; ...)

Deterministic loss-free multiwire summaries, `(mean_x, sigma_x, mean_y, sigma_y)`
in mm, one observation per row. `inputs` has the ordering from `bta_input_names()`.
This is the smooth function to compose with AD/OED. Its nested AD includes mixed
second tracking derivatives. Hard-aperture particle counting and noise are
provided separately by `simulate_bta`.
"""
function bta_measurements(
        inputs::AbstractVector;
        op=default_bta_operating_point(), observations=PROFILE_NAMES, moment_order=1,
        z0=0.0,
        quad_factors=ones(15), quad_dx=zeros(15), quad_dy=zeros(15),
        calibration=DEFAULT_BTA_CALIBRATION
    )
    length(inputs) == 20 + length(BEAM_PARAMETER_NAMES) ||
        throw(DimensionMismatch("wrong number of magnet and beam inputs"))
    state = UncoupledBeamState(inputs[21:end])
    beam = beam_moments_from_twiss(state; z0)
    setting = merge(op, (; quad_currents=inputs[1:15], corrector_currents=inputs[16:20]))
    predict_profiles(beam; op=setting, observations, moment_order, quad_factors,
                     quad_dx, quad_dy, calibration)
end

bta_measurements(quad_currents, corrector_currents, beam_state; kwargs...) =
    bta_measurements(pack_bta_inputs(quad_currents, corrector_currents, beam_state); kwargs...)

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

function _gaussian_particles(state::UncoupledBeamState, normals; z0)
    L = _twiss_factor(; _beam_optics(state)...)
    normals * transpose(L) .+ transpose(_canonical_mean(state, z0))
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
    simulate_bta(beam_state=UncoupledBeamState(); op=default_bta_operating_point(), ...)

Simulate an entrance `UncoupledBeamState` (or an equivalent keyword NamedTuple)
at operating point `op`, which supplies all magnet currents. The default
tracks 1000 Gaussian particles.

Returns clean/noisy loss-free moment observations and their 16 × 34 Jacobian,
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
function simulate_bta(beam_state::Union{UncoupledBeamState,NamedTuple}=UncoupledBeamState();
        op=default_bta_operating_point(), observations=PROFILE_NAMES,
        noise_std_mm=0.1, rng::AbstractRNG=Random.default_rng(), nparticles::Integer=1000,
        standard_particles=nothing, compute_jacobian::Bool=true, chunksize::Integer=4,
        moment_order=1, z0=0.0,
        quad_factors=ones(15), quad_dx=zeros(15), quad_dy=zeros(15),
        calibration=DEFAULT_BTA_CALIBRATION, apertures::Bool=true
    )
    nparticles > 0 || throw(ArgumentError("nparticles must be positive"))
    chunksize > 0 || throw(ArgumentError("chunksize must be positive"))
    state = beam_state isa NamedTuple ? UncoupledBeamState(beam_state) : beam_state
    inputs = Float64.(pack_bta_inputs(op.quad_currents, op.corrector_currents, state))
    sigma = _summary_noise_scale(noise_std_mm, length(observations))
    isfinite(z0) || throw(ArgumentError("longitudinal position must be finite"))
    if !isnothing(standard_particles)
        size(standard_particles) == (nparticles, 5) ||
            throw(DimensionMismatch("standard_particles must be ($nparticles, 5)"))
        all(isfinite, standard_particles) || throw(ArgumentError("standard particles must be finite"))
    end
    options = (; op, observations, moment_order, z0, quad_factors, quad_dx, quad_dy, calibration)
    clean = bta_measurements(inputs; options...)
    jacobian = compute_jacobian ?
        _measurement_value_jacobian(inputs, options, Val(min(chunksize, length(inputs)))) : nothing
    perturbation = sigma .* randn(rng, size(clean)...)
    noisy = clean + perturbation
    normals = isnothing(standard_particles) ? randn(rng, nparticles, 5) : Float64.(standard_particles)
    # The particle path is deliberately primal; differentiate bta_measurements.
    primal_state = UncoupledBeamState(inputs[21:end])
    initial_particles = _gaussian_particles(primal_state, normals; z0)
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

# Preserve positional current overrides for existing callers.
function simulate_bta(quad_currents, corrector_currents, beam_state=UncoupledBeamState();
                      op=default_bta_operating_point(), kwargs...)
    simulate_bta(beam_state; op=merge(op, (; quad_currents, corrector_currents)), kwargs...)
end
