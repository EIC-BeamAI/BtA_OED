mutable struct BtaWorkspace{L,E,I,O,C}
    lattice::L
    elements::E
    indices::I
    observations::O
    calibration::C
    apertures::Bool
    revision::UInt
end

"""Prepare one independently owned lattice for repeated evaluations."""
function prepare_bta(op=default_bta_operating_point(); observations=PROFILE_NAMES,
                     apertures=false, calibration=DEFAULT_BTA_CALIBRATION, settings...)
    names = Tuple(Symbol.(observations))
    isempty(names) && throw(ArgumentError("at least one observation is required"))
    lattice = bta_lattice(op; apertures, calibration, settings...)
    elements = Dict(Symbol(ele.name) => ele for ele in lattice.line)
    BtaWorkspace(lattice, elements, _observation_indices(lattice, names), names,
                 calibration, apertures, UInt(0))
end

"""Apply absolute settings to a workspace. Use one workspace per worker."""
function set_operating_point!(work::BtaWorkspace, op; quad_factors=ones(15),
                              quad_dx=zeros(15), quad_dy=zeros(15), calibration=work.calibration)
    op.species_ref == work.lattice.species_ref ||
        throw(ArgumentError("prepare a separate workspace for a different reference species"))
    length(quad_dx) == length(quad_dy) == length(QUAD_NAMES) ||
        throw(DimensionMismatch("one offset per quadrupole and plane is required"))
    work.revision += 1
    set_bta_fields!(work.lattice, op; quad_factors, calibration, work.elements, reuse_groups=true)
    set_bta_misalignments!(work.lattice; quad_dx, quad_dy, work.elements)
    work
end

function _prepared_tensors(work, mean, order)
    work.apertures && throw(ArgumentError("moment predictions require an aperture-free workspace"))
    _transport_tensors(work.lattice, mean, work.indices, work.observations; hessians=order == 2)
end

"""Return a callable beam predictor. Supply fixed_mean to cache transport for covariance-only scans."""
function prepared_prediction(work::BtaWorkspace; output=:profiles, moment_order=1,
                             fixed_mean=nothing)
    output in (:centroid, :profiles, :moments) || throw(ArgumentError("unknown prediction output"))
    moment_order in (1, 2) || throw(ArgumentError("moment_order must be 1 or 2"))
    work.apertures && throw(ArgumentError("moment predictions require an aperture-free workspace"))
    centroid_only = output == :centroid && moment_order == 1
    evaluate(beam, tensors) = if centroid_only
        reshape(_transport_coordinates(work.lattice, beam.mean, work.indices), 6, :)
    else
        moments = _gaussian_closure(tensors, beam.covariance, moment_order)
        output == :profiles ? _profile_summaries(moments) :
            output == :centroid ? moments.means : moments
    end
    if isnothing(fixed_mean)
        return beam -> evaluate(beam, centroid_only ? nothing : _prepared_tensors(work, beam.mean, moment_order))
    end
    revision = work.revision
    mean = copy(fixed_mean)
    cached = centroid_only ? evaluate(BeamMoments(mean, zeros(6, 6); check=false), nothing) :
        _prepared_tensors(work, mean, moment_order)
    beam -> begin
        work.revision == revision || throw(ArgumentError("lattice changed; prepare the cached predictor again"))
        beam.mean == mean && !any(x -> x isa ForwardDiff.Dual, beam.mean) ||
            throw(ArgumentError("cached prediction requires the fixed, inactive centroid"))
        centroid_only ? copy(cached) : evaluate(beam, cached)
    end
end

"""Prepare an in-place Jacobian for a scalar or array-valued function of the selected inputs."""
function prepare_differential(f, inputs::AbstractVector; chunksize=min(length(inputs), 4))
    1 <= chunksize <= length(inputs) || throw(ArgumentError("invalid Jacobian chunk size"))
    _prepare_differential(f, inputs, Val(chunksize))
end

_differential_output(value::Real) = [value]
_differential_output(value::AbstractArray) = vec(value)

function _prepare_differential(f, inputs, ::Val{C}) where {C}
    flatten = p -> _differential_output(f(p))
    config = ForwardDiff.JacobianConfig(flatten, inputs, ForwardDiff.Chunk{C}())
    rows = length(flatten(inputs))
    (rows=rows, columns=length(inputs),
     jacobian! = (jac, p) -> ForwardDiff.jacobian!(jac, flatten, p, config))
end

"""Track an N × 6 particle matrix through the prepared lattice."""
function track_particles(work::BtaWorkspace, coordinates::AbstractMatrix)
    size(coordinates, 2) == 6 && size(coordinates, 1) > 0 ||
        throw(DimensionMismatch("expected a nonempty N × 6 particle matrix"))
    bunch = Bunch(copy(coordinates); species=work.lattice.species_ref,
                  p_over_q_ref=work.lattice.p_over_q_ref)
    BeamTracking.track!(bunch, work.lattice; rf_on=false, use_explicit_SIMD=false,
                        use_cpu_multithreading=false)
    bunch
end

"""Return the surviving particle fraction using a prepared lattice."""
function transmission!(work::BtaWorkspace, coordinates::AbstractMatrix)
    bunch = track_particles(work, coordinates)
    count(==(BeamTracking.STATE_ALIVE), bunch.coords.state) / size(coordinates, 1)
end
