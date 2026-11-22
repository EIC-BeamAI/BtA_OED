# BeamlineExperiments plotting module

This directory contains a Julia module and examples, not a separate package.
All dependencies and source pins are in BeamlineExperiments's root Project/Manifest.

## Installation and interactive viewer

From the repository root:

```sh
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. -i plotting/examples/apertures.jl
```

The parser and bundled JLL are pinned to PALSParserJ's `feat/jll-packaging`
commit `59ac1765398c14723ce6fb6e1ddcdf34eafa5f4e`. The JLL is a **remote Git
subdirectory dependency** (`contrib/PALSParserCpp_jll`), not a local checkout.
Julia downloads the platform binary from permanent upstream URLs and verifies
its checksum. No native compilation, manual artifact installation, GitHub CLI,
or separate bootstrap is needed. Native targets: macOS ARM64/Intel, Linux
x86_64 glibc, and Windows x86_64.

In a root-project Julia REPL:

```julia
using GLMakie
include("plotting/BeamlineExperimentsPlotting.jl")
using .BeamlineExperimentsPlotting
GLMakie.activate!()
view = interactive_apertures()
# Or use a newly exported PALS file:
view = interactive_apertures("my-bta.yaml")
```

Scroll to zoom, right-drag to pan, left-drag to choose a zoom rectangle.
The 1–100× slider changes only the transverse aperture offsets; the curved
reference centerline stays fixed. **Reset view** fits the current envelope.
The lower panel shows the actual horizontal limits in millimetres.

Bands illustrate nominal local-x apertures along the elements' lengths;
body misalignments are omitted. Dots mark the aperture-check locations.
Missing, inactive, or unbounded apertures remain gaps. This is not a beam
size or transmission calculation.

The module does not select a backend on import. Use GLMakie on a desktop with
OpenGL, or CairoMakie for headless rendering. Start a fresh session after changing
parser versions; remove legacy PALS_PARSER_CPP_LIB/PALS_PARSER_CPP_DIR overrides
to let the JLL choose its binary. The root shell helper clears those overrides.

## Floor plan and static image

```julia
floor = floor_view()
display(floor)
```

```sh
julia --project=. plotting/examples/render_apertures.jl
julia --project=. plotting/examples/render_apertures.jl my-bta.yaml output.png
```

Output defaults to `plotting/output/bta-apertures.png`. The viewer uses
`src/bta_lattice/bta.pals.yaml` by default.

Run the headless plotting checks with `julia --project=. test/plotting.jl`.
For the optional desktop OpenGL check, run
`julia --project=. test/plotting_glmakie.jl`.
