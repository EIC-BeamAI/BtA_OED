module BtaUQ

using Beamlines
using BeamTracking
using ForwardDiff
using LinearAlgebra
using Random
using SpecialFunctions: erf
import ChainRulesCore

export QUAD_NAMES, QUAD_ELEMENT_NAMES, CORRECTOR_NAMES, PROFILE_NAMES,
    DEFAULT_BTA_CALIBRATION, default_bta_operating_point, bta_fields,
    bta_lattice, set_bta_fields!, set_bta_misalignments!,
    BeamMoments, default_beam_moments, beam_moments_from_twiss,
    transport_coordinates, transport_tensors, predict_moments, predict_profiles,
    profile_vector, profile_value_and_jacobian, profile_jacobians,
    profile_probabilities, track_particles, transmission,
    BatchedForwardSensitivity, PreparedBtaBatch, BtaBatchProblem,
    prepare_bta_batch, batch_problem, chain_problems, bta_loglikelihood,
    value_and_gradient!, factor_batch_loglikelihoods, factor_batch_value_and_gradient!,
    pack_parameters, unpack_parameters, parameter_names,
    UncoupledBeamState, uncoupled_beam_parameters, pack_bta_inputs,
    bta_input_names, bta_measurements, simulate_bta

include("bta_lattice/bta_conversions.jl")
include("bta_lattice/bta_lattice.jl")
include("bta_lattice/bta_setting.jl")
include("bta_moments.jl")
include("bta_simulation.jl")
#include("bta_uq_preparation.jl")
#include("bta_uq_inference.jl")

end
