module BeamlineExperiments

using Beamlines
using BeamTracking
using ForwardDiff
using LinearAlgebra
using Random
using SpecialFunctions: erf

export QUAD_NAMES, QUAD_ELEMENT_NAMES, CORRECTOR_NAMES, SEPTUM_NAMES, PROFILE_NAMES,
    DEFAULT_BTA_CALIBRATION, default_bta_operating_point, bta_fields,
    bta_lattice, set_bta_fields!, set_bta_misalignments!,
    BeamMoments, default_beam_moments, beam_moments_from_twiss,
    transport_coordinates, transport_tensors, predict_moments, predict_profiles,
    profile_vector, profile_value_and_jacobian, profile_jacobians,
    profile_probabilities, track_particles, transmission,
    prepare_bta, set_operating_point!, prepared_prediction, prepare_differential, transmission!,
    BEAM_PARAMETER_NAMES, beam_decoder, machine_decoder, prepare_bta_batch, bta_loglikelihood,
    AGS_INJECTION_TARGET, injection_mismatch, aggregate_beam, control_box,
    UncoupledBeamState, uncoupled_beam_parameters, pack_bta_inputs,
    bta_input_names, bta_measurements, simulate_bta

include("bta_lattice/bta_conversions.jl")
include("bta_lattice/bta_lattice.jl")
include("bta_lattice/bta_setting.jl")
include("transport/moments.jl")
include("transport/simulation.jl")
include("transport/workspace.jl")
include("inference/preparation.jl")
include("inference/likelihood.jl")
include("control/injection.jl")
include("control/boxes.jl")

end
