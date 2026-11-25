# Supply order follows the physical line. QH2A/B have one shared calibration.
const QUAD_NAMES = (:qv1, :qh2, :qv3, :qh4, :qv5, :qh6, :qv7, :qh8,
                    :qv9, :qh10, :qv11, :qh12, :qv13, :qh14, :qv15)
const QUAD_ELEMENT_NAMES = ((:qv1,), (:qh2a, :qh2b), (:qv3,), (:qh4,),
    (:qv5,), (:qh6,), (:qv7,), (:qh8,), (:qv9,), (:qh10,), (:qv11,),
    (:qh12,), (:qv13,), (:qh14,), (:qv15,))
const CORRECTOR_NAMES = (:dv007, :dv030, :dh127, :dh158, :dv181)
const SEPTUM_NAMES = (:f6, :l20)
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
    f6=(-1.0196e-4,), # AGS Studies Report 283,
    l20=(6.0149e-5,), # AP-473 gap estimate, AP-217 effective length
)

"""Currents matched to the nominal AGS injection target for the default entrance beam."""
function default_bta_operating_point(; kwargs...)
    defaults = (
        quad_currents=(253.17947, 405.17948, 535.65982, 555.94257, 150.49952,
                       839.84366, 727.11538, 73.97488, 490.55114, 445.84191,
                       291.59416, 74.43011, 234.96670, 273.88004, 146.92053),
        dipole_currents=(218.51883, 3393.67603, 100.56230, 419.88626),
        septum_currents=(9196.38385, 14598.37000), # F6, L20
        corrector_currents=(-0.68541, 0.31214, 5.53796, 19.81779, -0.28446),
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
The `septum_fields` entries are integrated dipole fields in T·m. Both magnets
in each pair use the same conductor current from `op.septum_currents`.
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
    currents = op.septum_currents
    length(currents) == 2 || throw(DimensionMismatch("two septum currents required (F6, L20)"))
    all(isfinite, currents) || throw(ArgumentError("septum currents must be finite"))
    f6 = only(calibration.f6) * currents[1]
    l20 = only(calibration.l20) * currents[2]
    septum_fields = (
        dhf6a=f6 * 49.213 / (49.213 + 52.243),
        dhf6b=f6 * 52.243 / (49.213 + 52.243),
        dhl20a=l20 * (44.2337 + 0.920) / 88.0037,
        dhl20b=l20 * (41.930 + 0.920) / 88.0037,
    )
    (; quadrupoles, dipoles, correctors, septum_fields)
end
