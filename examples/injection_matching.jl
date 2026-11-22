using BeamlineExperiments

box = control_box(:quadrupole)
beam = beam_moments_from_twiss(UncoupledBeamState(; centroid=box.centroid, sigma_delta=5e-4))
work = prepare_bta(box.op; observations=(:patch_l20_3,))
predict = prepared_prediction(work; output=:moments)

# Differentiate the one-shot injection objective with respect to two quad currents.
function score(currents)
    quads = [i <= 2 ? currents[i] : box.op.quad_currents[i] for i in 1:15]
    set_operating_point!(work, merge(box.op, (; quad_currents=quads)))
    result = predict(beam)
    terminal = BeamMoments(result.means[:, 1], result.covariances[:, :, 1]; check=false)
    injection_mismatch(terminal; epsilon_x=1e-6, epsilon_y=1e-6).reward
end

currents = collect(box.op.quad_currents[1:2])
diff = prepare_differential(score, currents)
jacobian = zeros(diff.rows, diff.columns)
diff.jacobian!(jacobian, currents)
@show score(currents) jacobian
