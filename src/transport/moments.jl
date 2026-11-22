"""Mean and covariance in canonical `(x, px, y, py, z, delta)` coordinates, SI units."""
struct BeamMoments{M,C}
    mean::M
    covariance::C
    function BeamMoments(mean::AbstractVector, covariance::AbstractMatrix; check=true)
        length(mean) == 6 || throw(DimensionMismatch("mean must have six coordinates"))
        size(covariance) == (6, 6) || throw(DimensionMismatch("covariance must be 6 × 6"))
        if check
            all(isfinite, mean) && all(isfinite, covariance) ||
                throw(ArgumentError("beam moments must be finite"))
            isapprox(covariance, transpose(covariance); atol=1e-15, rtol=1e-12) ||
                throw(ArgumentError("covariance must be symmetric"))
            primal = _primal.(covariance)
            minimum(eigvals(Symmetric(primal))) >= -1e-12 * max(opnorm(primal), eps()) ||
                throw(ArgumentError("covariance must be positive semidefinite"))
        end
        new{typeof(copy(mean)),typeof(copy(covariance))}(copy(mean), copy(covariance))
    end
end

_primal(x) = x
_primal(x::ForwardDiff.Dual) = _primal(ForwardDiff.value(x))

"""Build entrance moments from uncoupled betatron optics and linear dispersion.

The transverse centroid and momentum offset are independent inputs. Dispersion
adds correlations around that centroid through the RMS relative momentum spread.
Pass a full canonical mean with `mean` or use `centroid`, `delta0`, and `z0`.
"""
function beam_moments_from_twiss(; centroid=(0.0, 0.0, 0.0, 0.0), delta0=0.0, z0=0.0,
                                mean=[centroid..., z0, delta0], optics...)
    L = _twiss_factor(; optics...)
    BeamMoments(mean, L * transpose(L))
end

function _twiss_factor(; epsilon_x=1e-6, beta_x=12.6, alpha_x=1.8,
                         epsilon_y=1e-6, beta_y=4.1, alpha_y=-0.65,
                         eta_x=2.7, eta_px=-0.42, sigma_delta=1e-4)
    beta_x > 0 && beta_y > 0 && epsilon_x > 0 && epsilon_y > 0 && sigma_delta >= 0 ||
        throw(ArgumentError("beta/emittance must be positive and momentum spread nonnegative"))
    T = promote_type(typeof.((epsilon_x, beta_x, alpha_x, epsilon_y, beta_y, alpha_y,
        eta_x, eta_px, sigma_delta))...)
    L = zeros(T, 6, 5)
    L[1, 1] = sqrt(epsilon_x * beta_x)
    L[2, 1] = -sqrt(epsilon_x / beta_x) * alpha_x
    L[2, 2] = sqrt(epsilon_x / beta_x)
    L[3, 3] = sqrt(epsilon_y * beta_y)
    L[4, 3] = -sqrt(epsilon_y / beta_y) * alpha_y
    L[4, 4] = sqrt(epsilon_y / beta_y)
    L[:, 5] .= sigma_delta .* (eta_x, eta_px, 0, 0, 0, 1)
    L
end

"""
Illustrative entrance beam, not a measured beam estimate. The lattice starts
before the extraction-coordinate patch, where zero is not the beamline centre.
"""
default_beam_moments(; mean=[0.1013, 0.00174, 0.0, 0.0, 0.0, 0.0], kwargs...) =
    beam_moments_from_twiss(; mean, kwargs...)

function _observation_indices(lattice, names)
    symbols = Symbol.(names)
    length(unique(symbols)) == length(symbols) ||
        throw(ArgumentError("observation names must be unique"))
    by_name = Dict(Symbol(ele.name) => j for (j, ele) in enumerate(lattice.line))
    [name == :end ? length(lattice.line) : get(by_name, name, 0) > 0 ?
        by_name[name] : throw(ArgumentError("unknown observation $name")) for name in symbols]
end

function _active_mean(lattice, mean)
    # Promote coordinates to the active magnet/alignment/reference scalar type.
    # A zero-valued anchor carries outer AD tags without changing the trajectory.
    anchor = zero(mean[1])
    for ele in lattice.line
        bm = ele.BMultipoleParams
        if !isnothing(bm)
            for v in bm.n
                anchor += zero(v)
            end
        end
        ap = ele.AlignmentParams
        isnothing(ap) || (anchor += zero(ap.x_offset) + zero(ap.y_offset))
    end
    anchor += zero(lattice.p_over_q_ref)
    [x + anchor for x in mean]
end

function _transport_coordinates(lattice, mean, indices; promoted=false)
    coordinates = promoted ? collect(mean) : _active_mean(lattice, mean)
    state = reshape(coordinates, 1, 6)
    bunch = Bunch(state; species=lattice.species_ref, p_over_q_ref=lattice.p_over_q_ref)
    result = Vector{eltype(state)}(undef, 6length(indices))
    for j in 1:maximum(indices)
        element = lattice.line[j]
        BeamTracking.track!(bunch, element; rf_on=false, use_cpu_multithreading=false,
                            use_explicit_SIMD=false)
        only(bunch.coords.state) == BeamTracking.STATE_ALIVE ||
            error("reference trajectory failed at $(element.name); check centroid and settings")
        for k in eachindex(indices)
            if indices[k] == j
                copyto!(result, 6(k - 1) + 1, bunch.coords.v, 1, 6)
            end
        end
    end
    all(isfinite, result) || error("nonfinite trajectory; check the operating point")
    result
end

"""Unclipped coordinates at named observations, as a 6 × n matrix."""
function transport_coordinates(
        mean; 
        op=default_bta_operating_point(),
        quad_factors=ones(15), quad_dx=zeros(15), quad_dy=zeros(15),
        observations=PROFILE_NAMES, calibration=DEFAULT_BTA_CALIBRATION
    )
    length(mean) == 6 || throw(DimensionMismatch("six input coordinates required"))
    isempty(observations) && throw(ArgumentError("at least one observation required"))
    lattice = bta_lattice(op; apertures=false, quad_factors, quad_dx, quad_dy, calibration)
    indices = _observation_indices(lattice, observations)
    reshape(_transport_coordinates(lattice, mean, indices), 6, length(indices))
end

"""
Return coordinates, input-coordinate Jacobians and optionally coordinate Hessians.
Jacobians have shape `(6, 6, observations)`, Hessians `(6, 6, 6, observations)`.
All AD layers remain active, including mixed coordinate/magnet derivatives.
"""
function transport_tensors(mean; hessians=false, op=default_bta_operating_point(),
        quad_factors=ones(15), quad_dx=zeros(15), quad_dy=zeros(15),
        observations=PROFILE_NAMES, calibration=DEFAULT_BTA_CALIBRATION
    )
    isempty(observations) && throw(ArgumentError("at least one observation required"))
    lattice = bta_lattice(op; apertures=false, quad_factors, quad_dx, quad_dy, calibration)
    indices = _observation_indices(lattice, observations)
    _transport_tensors(lattice, mean, indices, observations; hessians)
end

function _transport_tensors(lattice, mean, indices, observations; hessians=false)
    f = x -> _transport_coordinates(lattice, x, indices; promoted=true)
    # Promote before seeding coordinate derivatives. Captured magnet duals must
    # be the value type of the inner dual, not a separately ordered peer tag.
    active_mean = _active_mean(lattice, mean)
    value = f(active_mean)
    J = ForwardDiff.jacobian(f, active_mean)
    n = length(indices)
    jacobians = permutedims(reshape(J, 6, n, 6), (1, 3, 2))
    H = if hessians
        dJ = ForwardDiff.jacobian(x -> vec(ForwardDiff.jacobian(f, x)), active_mean)
        permutedims(reshape(dJ, 6, n, 6, 6), (1, 3, 4, 2))
    else
        nothing
    end
    (; coordinates=reshape(value, 6, n), jacobians, hessians=H,
       observations=Tuple(Symbol.(observations)))
end

"""
    predict_moments(beam=default_beam_moments(); moment_order=1, ...)

Unclipped Gaussian moment closure. Order 1 uses `mu=f(mu0)`, `S=J*S0*J'`.
Order 2 adds `tr(H_i*S0)/2` to each mean and
`tr(H_i*S0*H_j*S0)/2` to the covariance (quadratic Taylor Gaussian closure).
The latter is not an exact nonlinear distribution propagation.

Differentiating order-1 moments with respect to magnets or the input centroid
requires second/mixed derivatives of tracking; nested AD retains these terms.
Order-2 moment gradients additionally require third derivatives.
"""
function predict_moments(beam::BeamMoments=default_beam_moments(); moment_order=1, kwargs...)
    moment_order in (1, 2) || throw(ArgumentError("moment_order must be 1 or 2"))
    tensors = transport_tensors(beam.mean; hessians=moment_order == 2, kwargs...)
    _gaussian_closure(tensors, beam.covariance, moment_order)
end

function _gaussian_closure(tensors, S, moment_order)
    means = [copy(tensors.coordinates[:, j]) for j in axes(tensors.coordinates, 2)]
    covariances = map(axes(tensors.coordinates, 2)) do k
        J = tensors.jacobians[:, :, k]
        C = J * S * transpose(J)
        if moment_order == 2
            HS = [tensors.hessians[i, :, :, k] * S for i in 1:6]
            means[k] = means[k] + [tr(h) / 2 for h in HS]
            C = C + [tr(HS[i] * HS[j]) / 2 for i in 1:6, j in 1:6]
        end
        (C + transpose(C)) / 2
    end
    (; means=reduce(hcat, means), covariances=cat(covariances...; dims=3),
       jacobians=tensors.jacobians, observations=tensors.observations, moment_order)
end

"""Rows are observations; columns `(mean_x, sigma_x, mean_y, sigma_y)`, in mm."""
function predict_profiles(beam::BeamMoments=default_beam_moments(); kwargs...)
    moments = predict_moments(beam; kwargs...)
    _profile_summaries(moments)
end

function _profile_summaries(moments)
    n = size(moments.means, 2)
    [column == 1 ? 1000 * moments.means[1, k] :
     column == 2 ? 1000 * sqrt(moments.covariances[1, 1, k]) :
     column == 3 ? 1000 * moments.means[3, k] :
                   1000 * sqrt(moments.covariances[3, 3, k])
     for k in 1:n, column in 1:4]
end

profile_vector(beam::BeamMoments=default_beam_moments(); kwargs...) =
    vec(predict_profiles(beam; kwargs...))

"""Jacobian with respect to an arbitrary caller-defined parameter vector/decoder."""
function profile_value_and_jacobian(parameters, decoder; kwargs...)
    f = p -> begin
        decoded = decoder(p)
        profile_vector(decoded.beam; decoded.settings..., kwargs...)
    end
    (; value=f(parameters), jacobian=ForwardDiff.jacobian(f, parameters))
end

"""Separate profile Jacobians for factors, centroid, covariance entries, and offsets."""
function profile_jacobians(
        beam::BeamMoments=default_beam_moments();
        quad_factors=ones(15), quad_dx=zeros(15), quad_dy=zeros(15), kwargs...
    )
    predict = (b, q, dx, dy) -> profile_vector(b; quad_factors=q, quad_dx=dx, quad_dy=dy, kwargs...)
    value = predict(beam, quad_factors, quad_dx, quad_dy)
    factors = ForwardDiff.jacobian(q -> predict(beam, q, quad_dx, quad_dy), quad_factors)
    centroid = ForwardDiff.jacobian(m ->
        predict(BeamMoments(m, beam.covariance; check=false), quad_factors, quad_dx, quad_dy), beam.mean)
    # Full vec(S) derivative: use symmetric perturbations so every AD direction
    # is a valid symmetric matrix variation, even for off-diagonal entries.
    covariance = ForwardDiff.jacobian(vec(beam.covariance)) do s
        S = reshape(s, 6, 6)
        predict(BeamMoments(beam.mean, (S + transpose(S)) / 2; check=false),
                quad_factors, quad_dx, quad_dy)
    end
    dx = ForwardDiff.jacobian(x -> predict(beam, quad_factors, x, quad_dy), quad_dx)
    dy = ForwardDiff.jacobian(y -> predict(beam, quad_factors, quad_dx, y), quad_dy)
    (; value, factors, centroid, covariance, quad_dx=dx, quad_dy=dy)
end

"""Gaussian wire-bin probabilities, excluding any charge/gain/background model."""
function profile_probabilities(edges, mean, sigma)
    sigma > 0 || throw(ArgumentError("sigma must be positive"))
    length(edges) >= 2 && issorted(edges) && all(diff(edges) .> 0) ||
        throw(ArgumentError("strictly increasing bin edges required"))
    cdf(x) = (1 + erf((x - mean) / (sqrt(2) * sigma))) / 2
    diff(cdf.(edges))
end

"""Hard-aperture particle tracking. Coordinates are N × 6."""
function track_particles(
        coordinates::AbstractMatrix; op=default_bta_operating_point(),
        quad_factors=ones(15), quad_dx=zeros(15), quad_dy=zeros(15),
        calibration=DEFAULT_BTA_CALIBRATION, apertures=true
    )
    size(coordinates, 2) == 6 || throw(DimensionMismatch("particles must be N × 6"))
    size(coordinates, 1) > 0 || throw(ArgumentError("at least one particle required"))
    lattice = bta_lattice(op; quad_factors, quad_dx, quad_dy, calibration, apertures)
    bunch = Bunch(copy(coordinates); species=lattice.species_ref, p_over_q_ref=lattice.p_over_q_ref)
    BeamTracking.track!(bunch, lattice; rf_on=false, use_explicit_SIMD=false,
                        use_cpu_multithreading=false)
    bunch
end

function transmission(coordinates::AbstractMatrix; kwargs...)
    bunch = track_particles(coordinates; kwargs...)
    count(==(BeamTracking.STATE_ALIVE), bunch.coords.state) / size(coordinates, 1)
end
