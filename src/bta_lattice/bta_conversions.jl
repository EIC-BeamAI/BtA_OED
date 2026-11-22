# Supply order follows the physical line. QH2A/B have one shared calibration.
const QUAD_NAMES = (:qv1, :qh2, :qv3, :qh4, :qv5, :qh6, :qv7, :qh8,
                    :qv9, :qh10, :qv11, :qh12, :qv13, :qh14, :qv15)
const QUAD_ELEMENT_NAMES = ((:qv1,), (:qh2a, :qh2b), (:qv3,), (:qh4,),
    (:qv5,), (:qh6,), (:qv7,), (:qh8,), (:qv9,), (:qh10,), (:qv11,),
    (:qh12,), (:qv13,), (:qh14,), (:qv15,))
const CORRECTOR_NAMES = (:dv007, :dv030, :dh127, :dh158, :dv181)
const PROFILE_NAMES = (:mw006, :mw060, :mw125, :mw166)

# Ascending powers, unnormalized integrated fields in SI units.
# Source decimal coefficients are represented in Float64 throughout.
const DEFAULT_BTA_CALIBRATION = (
    narrow=(4.363e-2, 9.266e-3, 6.309e-6, -8.609e-9, 2.226e-12),
    broad=(1.907e-2, 3.868e-3, 2.538e-6, -2.952e-9, 6.995e-13),
    short_dipole=(1.990e-3, 1.071e-3, 7.783e-7, -1.531e-9, 4.798e-13),
    strong_dipole_pair=(2.267e-3, 5.741e-4, 2.653e-8, -1.816e-11, 5.180e-15, -5.375e-19),
    long_dipole=(1.103e-2, 2.141e-3, 2.217e-6, -4.204e-9, 1.555e-12),
    dv007=(1.586e-4, 6.589e-4, -9.611e-8),
    corrector=(1.000e-4, 7.944e-4, -2.479e-6),
)

function default_bta_operating_point(; kwargs...)
    defaults = (
        quad_currents=(256.0, 406.0, 543.0, 556.0, 128.0, 797.0, 727.0, 0.0,
                       466.0, 460.0, 300.0, 72.0, 288.0, 0.0, 153.0),
        dipole_currents=(221.2, 3385.6, 106.6, 409.4), #(251.0, 3421.0, 112.0, 423.0),
        corrector_currents=(0.0, 0.0, 0.0, 0.0, 0.0),
        E_ref=4.5 * 938.272e6 / 1.79285,
        species_ref=Species("proton"),
    )
    unknown = setdiff(keys(kwargs), keys(defaults))
    isempty(unknown) || throw(ArgumentError("unknown operating-point keys: $unknown"))
    merge(defaults, (; kwargs...))
end

_horner(x, coefficients) = evalpoly(x, coefficients)

"""
    bta_fields(op; quad_factors=ones(15), calibration=DEFAULT_BTA_CALIBRATION)

Physical integrated fields, before rigidity normalization. As in Booster UQ,
each factor scales a complete calibrated quadrupole field, not just the linear
coefficient of its current polynomial. Both QH2 pieces receive the same factor.
"""
function bta_fields(op; quad_factors=ones(15), calibration=DEFAULT_BTA_CALIBRATION)
    length(op.quad_currents) == 15 || throw(DimensionMismatch("15 quad currents required"))
    length(op.dipole_currents) == 4 || throw(DimensionMismatch("4 dipole currents required"))
    length(op.corrector_currents) == 5 || throw(DimensionMismatch("5 corrector currents required"))
    length(quad_factors) == 15 || throw(DimensionMismatch("15 quad factors required"))
    quadrupoles = map(1:15) do j
        coefficients = j <= 4 || j == 15 ? calibration.narrow : calibration.broad
        sign = isodd(j) ? -1 : 1
        sign * quad_factors[j] * _horner(op.quad_currents[j], coefficients)
    end
    dipoles = (
        _horner(op.dipole_currents[1], calibration.short_dipole),
        _horner(op.dipole_currents[2], calibration.strong_dipole_pair),
        _horner(op.dipole_currents[3], calibration.short_dipole),
        -_horner(op.dipole_currents[4], calibration.long_dipole),
    )
    correctors = map(1:5) do j
        _horner(op.corrector_currents[j], j == 1 ? calibration.dv007 : calibration.corrector)
    end
    (; quadrupoles, dipoles, correctors)
end
