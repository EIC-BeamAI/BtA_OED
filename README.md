# BtA_OED

A Julia testbed for **uncertainty quantification, optimal experiment design, and control** on Brookhaven’s Booster-to-AGS beam transfer line. The Julia package is `BeamlineExperiments`.

Think of the line as a system with adjustable magnet currents, an imperfectly known input beam, and four measurement stations. Each station measures the beam’s center and width in both transverse directions. The project uses this setup to study three questions:

- **Inference.** What can we learn about the incoming beam and magnet calibration from these measurements?
- **Experiment design.** Which magnet settings would make the next measurements most informative?
- **Control.** Which settings deliver a beam closest to the injection target, and how much beam survives along the way?

The package provides differentiable measurement predictions, a configurable likelihood, particle tracking for estimating losses, and a downstream matching score. These are building blocks for experiment design and control algorithms.

## Quick start

From the repository root, run `julia --project=. -e 'using Pkg; Pkg.instantiate()'` with Julia 1.12. Then start Julia with `julia --project=.` and run:

```julia
using BeamlineExperiments, Random

box = control_box(:quadrupole)
op = box.op
beam = UncoupledBeamState(; centroid=box.centroid)

r = simulate_bta(beam; op, rng=MersenneTwister(42))
score = injection_mismatch(beam; op)
@show size(r.jacobian) r.fractional_beam_loss score.reward
```

The Jacobian describes how the predicted measurements change with magnet currents and entrance beam parameters in the loss-free model. Particle tracking estimates the fraction lost to physical apertures separately.

See [profile inference](examples/profile_inference.jl) for a synthetic-data likelihood and [injection matching](examples/injection_matching.jl) for a differentiable control objective.
