"""Uncoupled AGS target in canonical coordinates at the exit of patch_l20_3."""
const AGS_INJECTION_TARGET = (
    centroid=(0.05, 0.024, 0.0, 0.0),
    x=(beta=19.841, alpha=-1.5733,
       dispersion=(-2.066, -0.1474)),
    y=(beta=11.727, alpha=1.026,
       dispersion=(0.0, 0.0)),
)

"""Pool physical shots, optionally weighting by surviving charge."""
function aggregate_beam(shots; weights=ones(length(shots)))
    length(weights) == length(shots) || throw(DimensionMismatch("one weight per shot required"))
    all(w -> isfinite(w) && w >= 0, weights) && sum(weights) > 0 ||
        throw(ArgumentError("shot weights must be finite and nonnegative with positive sum"))
    indices = findall(>(0), weights)
    length(indices) == 1 && return shots[only(indices)]
    weight = weights ./ sum(weights)
    mean = sum(weight[i] * shots[i].mean for i in indices)
    covariance = sum(indices) do i
        offset = shots[i].mean - mean
        weight[i] * (shots[i].covariance + offset * transpose(offset))
    end
    BeamMoments(mean, covariance; check=false)
end

"""Phase-mixed excess squared size and reward at AGS injection.

Pass a single deterministic beam or pool physical shots with aggregate_beam.
Nominal geometric emittances set the reference sizes. For equal actual and
nominal betatron emittances, the envelope term is the amplitude term of Eq. (3).
Using the actual covariance also captures emittance changes in pooled beams.
"""
function injection_mismatch(beam::BeamMoments; epsilon_x, epsilon_y,
                            target=AGS_INJECTION_TARGET, sigma_pz=5e-4, horizontal_dispersion=nothing)
    all(x -> isfinite(x) && x > 0, (epsilon_x, epsilon_y)) ||
        throw(ArgumentError("nominal emittances must be finite and positive"))
    isfinite(sigma_pz) && sigma_pz >= 0 || throw(ArgumentError("invalid momentum spread"))
    function plane(indices, optics, epsilon, horizontal)
        β, α = optics.beta, optics.alpha
        γ = (1 + α^2) / β
        action(v) = γ*v[1]^2 + 2α*v[1]*v[2] + β*v[2]^2
        C = beam.covariance
        variance = C[6, 6]
        d = C[indices, 6]
        if horizontal && variance == 0 && sigma_pz > 0 && isnothing(horizontal_dispersion)
            throw(ArgumentError("supply horizontal_dispersion when beam momentum variance is zero"))
        end
        dispersion = !horizontal ? zeros(2) :
            !isnothing(horizontal_dispersion) ? collect(horizontal_dispersion) :
            variance > 0 ? d / variance : zeros(2)
        B = C[indices, indices]
        if horizontal
            B = B - d * transpose(dispersion) - dispersion * transpose(d) +
                variance * dispersion * transpose(dispersion)
        end
        delta_dispersion = dispersion - collect(optics.dispersion)
        offset = beam.mean[indices] - collect(target.centroid)[indices] -
                 collect(optics.dispersion) * beam.mean[6]
        envelope = (γ*B[1, 1] + 2α*B[1, 2] + β*B[2, 2]) / (2epsilon) - 1
        dispersive = sigma_pz^2 * action(delta_dispersion) / (2epsilon)
        centroid = action(offset) / (2epsilon)
        (; envelope, dispersion=dispersive, centroid,
           excess=envelope + dispersive + centroid, nominal_size2=β*epsilon)
    end
    x = plane(1:2, target.x, epsilon_x, true)
    y = plane(3:4, target.y, epsilon_y, false)
    cost = (x.nominal_size2*x.excess + y.nominal_size2*y.excess) /
           (x.nominal_size2 + y.nominal_size2)
    (; x, y, cost, reward=-cost)
end

"""Track one entrance beam to the injection plane and score its loss-free moments.

The entrance emittances supply the nominal reference by default. Override them
with `epsilon_x` and `epsilon_y` to compare against a fixed external reference.
For repeated evaluations, use a prepared predictor and score the terminal moments.
"""
function injection_mismatch(entrance::UncoupledBeamState; op=default_bta_operating_point(),
                            epsilon_x=entrance.epsilon_x, epsilon_y=entrance.epsilon_y,
                            z0=0.0, moment_order=1, kwargs...)
    beam = beam_moments_from_twiss(entrance; z0)
    m = predict_moments(beam; op, observations=(:patch_l20_3,), moment_order)
    terminal = BeamMoments(m.means[:, 1], m.covariances[:, :, 1]; check=false)
    injection_mismatch(terminal; epsilon_x, epsilon_y, kwargs...)
end
