using Test
using LinearAlgebra
using ForwardDiff
using Beamlines
using BeamTracking
using BeamlineExperiments

setentry(x::Tuple, value, index) = Base.setindex(x, value, index)
setentry(x::AbstractArray, value, index) = (y = x .+ zero(value); y[index] = value; y)

central(f, x, h) = (f(x + h) - f(x - h)) / (2h)

@testset "BtA supply wiring and injection endpoint" begin
    op = default_bta_operating_point()
    lattice = bta_lattice(op)
    @test last(lattice.line).name == "patch_l20_3"
    beam = default_beam_moments()
    @test transport_coordinates(beam.mean; observations=(:end,)) ==
          transport_coordinates(beam.mean; observations=(:patch_l20_3,))
    # Previously ignored supply settings must now reach every corresponding piece.
    changed = default_bta_operating_point(
        dipole_currents=(260.0, 3421.0, 120.0, 423.0),
        septum_scales=(0.98, 1.02),
        quad_currents=setentry(op.quad_currents, 10.0, 14))
    moved = bta_lattice(changed)
    byname(bl) = Dict(Symbol(e.name) => e for e in bl.line)
    a, b = byname(lattice), byname(moved)
    septa = (:dhf6a, :dhf6b, :dhl20a, :dhl20b)
    @test [b[name].Kn0L for name in septa] ≈
          [a[name].Kn0L for name in septa] .* [0.98, 0.98, 1.02, 1.02]
    for name in (:dh1_u, :dh1_d, :dh4_u, :dh4_d, :qh14)
        @test a[name].BMultipoleParams.n != b[name].BMultipoleParams.n
    end
    factors = ones(15)
    factors[2] = 1.02
    scaled = byname(bta_lattice(op; quad_factors=factors))
    @test scaled[:qh2a].Bn1L ≈ 1.02a[:qh2a].Bn1L
    @test scaled[:qh2b].Bn1L ≈ 1.02a[:qh2b].Bn1L
    shifts = fill(1e-4, 15)
    set_bta_misalignments!(lattice; quad_dx=shifts)
    set_bta_misalignments!(lattice; quad_dx=shifts)
    @test byname(lattice)[:qv1].x_offset ≈ 1e-4
    @test all(isnothing(e.ApertureParams) for e in bta_lattice(; apertures=false).line)
end

println("Running: Independent linear transport limits"); flush(stdout)
@testset "Independent linear transport limits" begin
    drift = Beamline([Drift(name="d", L=2.0)]; species_ref=Species("proton"), E_ref=3e9)
    f = m -> BeamlineExperiments._transport_coordinates(drift, m, [length(drift.line)])
    J = ForwardDiff.jacobian(f, zeros(6))
    expected = [1.0 2.0 0.0 0.0; 0.0 1.0 0.0 0.0;
                0.0 0.0 1.0 2.0; 0.0 0.0 0.0 1.0]
    @test J[1:4, 1:4] ≈ expected
    @test norm(J[1:4, 1:4] - Matrix{Float64}(I, 4, 4)) > 1 # negative control
    k, L = 0.7, 0.4
    quad = Beamline([Quadrupole(name="q", L=L, Kn1=k)];
                    species_ref=Species("proton"), E_ref=3e9)
    g = m -> BeamlineExperiments._transport_coordinates(quad, m, [length(quad.line)])
    Q = ForwardDiff.jacobian(g, zeros(6))
    w = sqrt(k)
    focusing = [cos(w*L) sin(w*L)/w; -w*sin(w*L) cos(w*L)]
    defocusing = [cosh(w*L) sinh(w*L)/w; w*sinh(w*L) cosh(w*L)]
    @test Q[1:2, 1:2] ≈ focusing
    @test Q[3:4, 3:4] ≈ defocusing
end

println("Running: Moments and mixed second derivatives"); flush(stdout)
@testset "Moments and mixed second derivatives" begin
    beam = default_beam_moments()
    @test minimum(eigvals(Symmetric(beam.covariance))) >= -1e-20
    profiles = predict_profiles(beam)
    @test size(profiles) == (4, 4)
    @test all(isfinite, profiles)
    @test all(profiles[:, [2, 4]] .> 0)
    moments = predict_moments(beam)
    @test all(isapprox(moments.covariances[:, :, j],
        moments.jacobians[:, :, j] * beam.covariance * moments.jacobians[:, :, j]') for j in 1:4)

    # Only changing covariance must change widths, while leaving order-1 centroids fixed.
    doubled = predict_profiles(BeamMoments(beam.mean, 2beam.covariance))
    @test doubled[:, [1, 3]] ≈ profiles[:, [1, 3]]
    @test doubled[:, [2, 4]] ≈ sqrt(2) .* profiles[:, [2, 4]]

    # The last two quadrupoles are downstream of every harp, but visible at :end.
    q = ones(15)
    q[14:15] .= 1.1
    @test predict_profiles(beam; quad_factors=q) == profiles
    @test !isapprox(predict_profiles(beam; quad_factors=q, observations=(:end,)),
                   predict_profiles(beam; observations=(:end,)))

    jac = profile_jacobians(beam; observations=(:mw060,))
    @test size(jac.factors) == (4, 15)
    @test size(jac.centroid) == (4, 6)
    @test size(jac.covariance) == (4, 36)
    for h in (1e-4, 3e-5)
        f = t -> profile_vector(beam; observations=(:mw060,),
                               quad_factors=setentry(ones(15), t, 1))
        fd = central(f, 1.0, h)
        @test isapprox(jac.factors[:, 1], fd; rtol=2e-5, atol=2e-5)
    end
    for (key, h) in ((:centroid, 1e-6), (:quad_dx, 1e-6), (:quad_dy, 1e-6))
        f = if key == :centroid
            t -> profile_vector(BeamMoments(setentry(beam.mean, beam.mean[1] + t, 1), beam.covariance);
                                 observations=(:mw060,))
        else
            t -> profile_vector(beam; observations=(:mw060,),
                    (key => setentry(zeros(15), t, 1),)...)
        end
        @test isapprox(getproperty(jac, key)[:, 1], central(f, 0.0, h);
                       rtol=3e-5, atol=2e-4)
    end
    # A frozen transport Jacobian would give zero for these width derivatives.
    @test norm(jac.factors[[2, 4], 1]) > 1e-3
    @test norm(jac.centroid[[2, 4], 1]) > 1e-3

    # A design current gradient traverses both polynomial conversion and tracking J.
    op = default_bta_operating_point()
    f = current -> profile_vector(beam; observations=(:mw060,),
        op=default_bta_operating_point(quad_currents=setentry(op.quad_currents, current, 1)))
    ad = ForwardDiff.derivative(f, op.quad_currents[1])
    @test isapprox(ad, central(f, op.quad_currents[1], 1e-2); rtol=2e-5, atol=2e-5)

    # Symmetric covariance direction, checked independently by a finite perturbation.
    direction = zeros(6, 6)
    direction[1, 2] = direction[2, 1] = 1.0
    f = t -> profile_vector(BeamMoments(beam.mean, beam.covariance + t * direction);
                             observations=(:mw060,))
    @test isapprox(jac.covariance * vec(direction), central(f, 0.0, 1e-10);
                   rtol=1e-5, atol=1e-2)
    @test sum(profile_probabilities([-Inf, 0.0, Inf], 0.0, 1.0)) ≈ 1
    @test profile_probabilities([-Inf, 0.0, Inf], 0.0, 1.0) ≈ [0.5, 0.5]
end

println("Running: Quadratic Gaussian closure and higher tracking derivatives"); flush(stdout)
@testset "Quadratic Gaussian closure and higher tracking derivatives" begin
    # Analytic scalar map y=x^2: E[y]=mu^2+s^2, Var[y]=4mu^2*s^2+2s^4.
    mu, variance = 0.3, 0.04
    coordinates = zeros(6, 1)
    coordinates[1] = mu^2
    jacobians = zeros(6, 6, 1)
    jacobians[1, 1, 1] = 2mu
    hessians = zeros(6, 6, 6, 1)
    hessians[1, 1, 1, 1] = 2
    S = zeros(6, 6)
    S[1, 1] = variance
    tensors = (; coordinates, jacobians, hessians, observations=(:analytic,))
    closed = BeamlineExperiments._gaussian_closure(tensors, S, 2)
    @test closed.means[1] ≈ mu^2 + variance
    @test closed.covariances[1, 1, 1] ≈ 4mu^2 * variance + 2variance^2
    linear = BeamlineExperiments._gaussian_closure(tensors, S, 1)
    @test abs(linear.covariances[1, 1, 1] - closed.covariances[1, 1, 1]) > 1e-3

    beam = default_beam_moments()
    tensors = transport_tensors(beam.mean; hessians=true, observations=(:mw006,))
    @test size(tensors.hessians) == (6, 6, 6, 1)
    @test isapprox(tensors.hessians[:, :, :, 1],
                   permutedims(tensors.hessians[:, :, :, 1], (1, 3, 2)); atol=1e-10)
    # Recover Hessian coefficients by a polynomial fit, at two fit ranges.
    for radius in (1e-3, 5e-4)
        nodes = collect(range(-radius, radius; length=11))
        samples = [transport_coordinates(setentry(beam.mean, beam.mean[2] + t, 2);
                    observations=(:mw006,))[1] for t in nodes]
        design = [t^power for t in nodes, power in 0:4]
        coefficients = design \ samples
        @test isapprox(2coefficients[3], tensors.hessians[1, 2, 2, 1];
                       rtol=1e-3, atol=1e-6)
    end
    second = predict_moments(beam; moment_order=2, observations=(:mw006,))
    @test all(isfinite, second.means)
    @test minimum(eigvals(Symmetric(second.covariances[:, :, 1]))) > -1e-15
    # Differentiating quadratic closure exercises a third tracking derivative.
    f = t -> profile_vector(BeamMoments(setentry(beam.mean, t, 2), beam.covariance);
                            moment_order=2, observations=(:mw006,))
    ad = ForwardDiff.derivative(f, beam.mean[2])
    @test isapprox(ad, central(f, beam.mean[2], 1e-5); rtol=2e-5, atol=2e-5)
end

@testset "Prepared inference and injection matching" begin
    beam = default_beam_moments()
    op = default_bta_operating_point()
    observations = (:mw060,)
    work = prepare_bta(op; observations)
    predict = prepared_prediction(work)
    data = reshape(predict(beam), 1, 1, 4) .+ 0.03
    likelihood = bta_loglikelihood(prepare_bta_batch([op]; observations), data,
        machine_decoder(:quad_factors, (1,)); noise_std_mm=0.1)
    f = t -> likelihood([t])
    @test f(1.0) ≈ -2 * (0.3^2 + log(2π) + 2log(0.1))
    @test ForwardDiff.derivative(f, 1.0) ≈ central(f, 1.0, 1e-6) rtol=5e-4 atol=1e-3
    cached = prepared_prediction(work; fixed_mean=beam.mean)
    changed = BeamMoments(beam.mean, 2beam.covariance)
    @test cached(changed) ≈ predict(changed)
    set_operating_point!(work, op; quad_dx=[0.001; zeros(14)])
    set_operating_point!(work, op)
    @test predict(beam) ≈ predict_profiles(beam; observations)
    @test_throws ArgumentError cached(changed)

    t = AGS_INJECTION_TARGET
    matched = beam_moments_from_twiss(; beta_x=t.x.beta, alpha_x=t.x.alpha,
        beta_y=t.y.beta, alpha_y=t.y.alpha, eta_x=t.x.dispersion[1],
        eta_px=t.x.dispersion[2], sigma_delta=5e-4)
    score(b) = injection_mismatch(b; epsilon_x=1e-6, epsilon_y=1e-6)
    @test score(matched).cost ≈ 0 atol=1e-14
    shots = [BeamMoments([0., p, 0., 0., 0., 0.], matched.covariance) for p in (-1e-4, 1e-4)]
    @test score(aggregate_beam(shots)).x.excess ≈ t.x.beta * 1e-8 / (2e-6)
    @test score(first(shots)).x.centroid ≈ t.x.beta * 1e-8 / (2e-6)
end

println("Running: Hard aperture tracking is a separate physical path"); flush(stdout)
@testset "Hard aperture tracking is a separate physical path" begin
    coordinates = zeros(3, 6)
    coordinates[:, 1] .= 10.0
    @test transmission(coordinates) == 0.0
    coordinates = repeat(permutedims(default_beam_moments().mean), 2)
    @test size(track_particles(coordinates; apertures=false).coords.v) == (2, 6)
end

include("scibmad_integration.jl")
include("simulation_interface.jl")

include("plotting.jl")
