function _normal_field(value, order)
    z = zero(value)
    BMultipoleParams([value], [z], [z], [order], [false], [true])
end

# Write a field again whenever either value carries an AD tag. Primal equality
# alone does not establish that two dual values have the same derivatives.
_same_primal(a, b) = !(a isa ForwardDiff.Dual || b isa ForwardDiff.Dual) && isequal(a, b)

function _set_normal_field!(element, value, order)
    group = element.BMultipoleParams
    if isnothing(group) || length(group.n) != 1 || !_same_primal(only(group.n), value)
        element.BMultipoleParams = _normal_field(value, order)
    end
end

"""Set complete multipole groups, preserving ForwardDiff's nested dual types."""
function set_bta_fields!(lattice::Beamline, op;
                         quad_factors=ones(15), calibration=DEFAULT_BTA_CALIBRATION,
                         elements=Dict(Symbol(ele.name) => ele for ele in lattice.line),
                         reuse_groups=false)
    fields = bta_fields(op; quad_factors, calibration)
    _same_primal(lattice.E_ref, op.E_ref) || (lattice.E_ref = op.E_ref)
    for (j, names) in enumerate(QUAD_ELEMENT_NAMES), name in names
        if reuse_groups
            _set_normal_field!(elements[name], fields.quadrupoles[j], 2)
        else
            elements[name].BMultipoleParams = _normal_field(fields.quadrupoles[j], 2)
        end
    end
    dipole_pieces = (
        (:dh1_u, 1, 17.125 / 33.375), (:dh1_d, 1, 16.25 / 33.375),
        (:dh2, 2, 1.0), (:dh3, 2, 1.0),
        (:dh4_u, 3, 0.5), (:dh4_d, 3, 0.5),
        (:dh5_u, 4, 0.358 / 1.134), (:dh5_m, 4, 0.416 / 1.134),
        (:dh5_d, 4, 0.360 / 1.134),
    )
    for (name, j, fraction) in dipole_pieces
        value = fields.dipoles[j] * fraction
        if reuse_groups
            _set_normal_field!(elements[name], value, 1)
        else
            elements[name].BMultipoleParams = _normal_field(value, 1)
        end
    end
    for (j, name) in enumerate(CORRECTOR_NAMES)
        field = fields.correctors[j]
        group = elements[name].BMultipoleParams
        if reuse_groups && !isnothing(group) && length(group.n) == 1 &&
                _same_primal(only(group.n), field)
            continue
        end
        z = zero(field)
        # Preserve the source's magnetic orientation, including its sign convention.
        tilt = j in (1, 2, 5) ? z + π / 2 : z
        elements[name].BMultipoleParams =
            BMultipoleParams([field], [z], [tilt], [1], [false], [true])
    end
    lattice
end

"""
Set supply-wise quad offsets relative to the surveyed source offsets.
The two QH2 pieces share one offset. Repeated calls do not accumulate shifts.
"""
function set_bta_misalignments!(lattice; quad_dx=zeros(15), quad_dy=zeros(15),
                                elements=Dict(Symbol(ele.name) => ele for ele in lattice.line))
    length(quad_dx) == length(quad_dy) == 15 ||
        throw(DimensionMismatch("15 offsets per plane required"))
    for (j, names) in enumerate(QUAD_ELEMENT_NAMES), name in names
        ele = elements[name]
        dx, dy = quad_dx[j], quad_dy[j]
        if dx isa ForwardDiff.Dual || dy isa ForwardDiff.Dual || !iszero(dx) || !iszero(dy)
            ele.AlignmentParams = AlignmentParams(x_offset=dx, y_offset=dy)
        elseif !isnothing(ele.AlignmentParams)
            ele.AlignmentParams = nothing
        end
    end
    lattice
end

"""
    bta_lattice(op=default_bta_operating_point(); apertures=true, ...)

Fresh, independently owned lattice. Hard apertures are retained for particle
tracking. Moment transport explicitly uses `apertures=false`: its moments are
unclipped, and are not moments of the surviving distribution.
"""
function bta_lattice(op=default_bta_operating_point(); apertures=true,
                     quad_factors=ones(15), quad_dx=zeros(15), quad_dy=zeros(15),
                     calibration=DEFAULT_BTA_CALIBRATION)
    lattice = Beamline(_bta_geometry(); E_ref=op.E_ref, species_ref=op.species_ref)
    set_bta_fields!(lattice, op; quad_factors, calibration)
    set_bta_misalignments!(lattice; quad_dx, quad_dy)
    if !apertures
        for ele in lattice.line
            ele.ApertureParams = nothing
        end
    end
    lattice
end
