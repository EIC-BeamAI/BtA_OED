using BeamlineExperiments, Random

# Estimate centroid survival versus scale around a box midpoint. This checks
# sampled currents, so it does not certify an entire box.
kind = isempty(ARGS) ? :current : Symbol(first(ARGS))
box = control_box(kind)
work = prepare_bta(box.op; observations=(:end,), apertures=true)
entrance = beam_moments_from_twiss(UncoupledBeamState(; centroid=box.centroid))
particles = reshape(entrance.mean, 1, 6)
lower, upper = collect(box.lower), collect(box.upper)
center, halfwidth = (lower + upper) / 2, (upper - lower) / 2
samples = 2rand(MersenneTwister(42), length(center), 256) .- 1

function survives(currents)
    op = kind == :current ?
        merge(box.op, (; dipole_currents=currents[1:4], corrector_currents=currents[5:9])) :
        merge(box.op, (; quad_currents=currents))
    set_operating_point!(work, op)
    transmission!(work, particles)
end

target = 0.20
scales = 0.125:0.125:1.5
fractions = map(scales) do scale
    sum(v -> survives(center + scale .* halfwidth .* v), eachcol(samples)) / size(samples, 2)
end
best = argmin(abs.(fractions .- target))
scale = scales[best]
@show scale fractions[best]
@show center - scale * halfwidth center + scale * halfwidth
