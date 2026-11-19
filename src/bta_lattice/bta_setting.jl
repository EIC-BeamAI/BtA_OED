function _normal_field(value, order)
    z = zero(value)
    BMultipoleParams([value], [z], [z], [order], [false], [true])
end

"""Set complete multipole groups, preserving ForwardDiff's nested dual types."""
function set_bta_fields!(lattice::Beamline, op;
                         quad_factors=ones(15), calibration=DEFAULT_BTA_CALIBRATION)
    lattice.E_ref = op.E_ref
    elements = Dict(Symbol(ele.name) => ele for ele in lattice.line)
    fields = bta_fields(op; quad_factors, calibration)
    for (j, names) in enumerate(QUAD_ELEMENT_NAMES), name in names
        elements[name].BMultipoleParams = _normal_field(fields.quadrupoles[j], 2)
    end
    dipole_pieces = (
        (:dh1_u, 1, 17.125 / 33.375), (:dh1_d, 1, 16.25 / 33.375),
        (:dh2, 2, 1.0), (:dh3, 2, 1.0),
        (:dh4_u, 3, 0.5), (:dh4_d, 3, 0.5),
        (:dh5_u, 4, 0.358 / 1.134), (:dh5_m, 4, 0.416 / 1.134),
        (:dh5_d, 4, 0.360 / 1.134),
    )
    for (name, j, fraction) in dipole_pieces
        elements[name].BMultipoleParams = _normal_field(fields.dipoles[j] * fraction, 1)
    end
    for (j, name) in enumerate(CORRECTOR_NAMES)
        field = fields.correctors[j]
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
function set_bta_misalignments!(lattice; quad_dx=zeros(15), quad_dy=zeros(15))
    length(quad_dx) == length(quad_dy) == 15 ||
        throw(DimensionMismatch("15 offsets per plane required"))
    elements = Dict(Symbol(ele.name) => ele for ele in lattice.line)
    for (j, names) in enumerate(QUAD_ELEMENT_NAMES), name in names
        ele = elements[name]
        if !(isapprox(quad_dx[j], zero(quad_dx[j])) && isapprox(quad_dy[j], zero(quad_dy[j])))
            ele.AlignmentParams = AlignmentParams(
                x_offset=quad_dx[j], y_offset=quad_dy[j]
            )
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
