"""Archived 20% centroid-survival scan ranges in amperes and their fixed settings."""
function control_box(kind::Symbol)
    centroid = (0.10132, 0.00174, 0.0, 0.0, 0.0)
    # Keep the scan baseline independent of the injection-matched default.
    baseline = default_bta_operating_point(;
        quad_currents=(256.0, 406.0, 543.0, 556.0, 128.0, 797.0, 727.0, 0.0,
                       466.0, 460.0, 300.0, 72.0, 288.0, 0.0, 153.0),
        dipole_currents=(221.2, 3385.6, 106.6, 409.4),
        septum_currents=(9262.65, 14380.0),
        corrector_currents=(0.0, 0.0, 0.0, 0.0, 0.0))
    if kind == :current
        return (; centroid, op=baseline,
            names=(:dh1, :dh2_dh3, :dh4, :dh5, CORRECTOR_NAMES...),
            lower=(208.26,3371.79,94.12,398.93,-12.23,-6.73,-13.56,-25.10,-15.57),
            upper=(227.89,3396.78,125.03,420.00,7.44,14.06,9.97,20.73,14.09))
    elseif kind == :quadrupole
        op = merge(baseline, (;
            dipole_currents=(221.21044953765244, 3385.62472653646, 106.62404879417201, 409.4269087684229)))
        return (; centroid, op, names=QUAD_NAMES,
            lower=(185.50,294.19,393.46,402.88,92.75,577.50,526.78,0.00,337.66,333.31,217.38,52.17,208.68,0.00,110.86),
            upper=(326.50,517.81,692.54,709.12,163.25,1016.50,927.22,440.64,594.34,586.69,382.62,91.83,367.32,440.64,195.14))
    end
    throw(ArgumentError("expected :current or :quadrupole"))
end
