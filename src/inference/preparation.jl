"""Decode selected entrance parameters by name or index in BEAM_PARAMETER_NAMES.

Unselected values stay at the supplied defaults. The longitudinal position z0
is fixed. Use a custom callable for transformed parameters.
"""
function beam_decoder(parameters; beam=UncoupledBeamState(), z0=0.0)
    selected = map(collect(parameters)) do parameter
        index = parameter isa Symbol ? findfirst(==(parameter), BEAM_PARAMETER_NAMES) : parameter
        index isa Integer && 1 <= index <= length(BEAM_PARAMETER_NAMES) ||
            throw(ArgumentError("unknown beam parameter $parameter"))
        Int(index)
    end
    allunique(selected) || throw(ArgumentError("duplicate beam parameters"))
    base = uncoupled_beam_parameters(beam)
    fixed_centroid = all(>(5), selected)
    nominal_mean = _canonical_mean(beam, z0)
    θ -> begin
        length(θ) == length(selected) || throw(DimensionMismatch("wrong beam parameter count"))
        values = similar(base, promote_type(eltype(base), eltype(θ)))
        values .= base
        values[selected] .= θ
        state = UncoupledBeamState(values)
        mean = fixed_centroid ? nominal_mean : _canonical_mean(state, z0)
        (; beam=beam_moments_from_twiss(state; mean))
    end
end

"""Build a decoder for selected quad factors or offsets, with a fixed entrance beam."""
function machine_decoder(field, indices; beam=default_beam_moments())
    field in (:quad_factors, :quad_dx, :quad_dy) ||
        throw(ArgumentError("expected quad_factors, quad_dx or quad_dy"))
    selected = collect(Int, indices)
    all(i -> 1 <= i <= length(QUAD_NAMES), selected) && allunique(selected) ||
        throw(ArgumentError("invalid quadrupole indices"))
    θ -> begin
        length(θ) == length(selected) || throw(DimensionMismatch("wrong machine parameter count"))
        values = fill(field == :quad_factors ? one(eltype(θ)) : zero(eltype(θ)), length(QUAD_NAMES))
        values[selected] .= θ
        (; beam, machine=NamedTuple{(field,)}((values,)))
    end
end

"""Prepare independent experiment predictors. The likelihood accepts any callable
returning (; beam) or (; beam, machine), where machine contains quad settings.
Use fixed_mean only for beam covariance inference with fixed machine settings.
"""
function prepare_bta_batch(operating_points; observations=PROFILE_NAMES,
                           fixed_mean=nothing, moment_order=1, calibration=DEFAULT_BTA_CALIBRATION)
    experiments = map(collect(operating_points)) do op
        work = prepare_bta(op; observations, calibration)
        (; op, work, predict=prepared_prediction(work; moment_order, fixed_mean))
    end
    isempty(experiments) && throw(ArgumentError("at least one experiment required"))
    (; experiments, fixed_mean)
end
