"""Return a conditional Gaussian log likelihood on a caller-defined latent vector.

Data shape is (experiments, monitors, 4). Standard deviations are in mm and can
be scalar or data-shaped. A mask can exclude unavailable measurements.
"""
function bta_loglikelihood(prepared, data, decode; noise_std_mm, mask=trues(size(data)))
    experiments = prepared.experiments
    shape = (length(experiments), length(first(experiments).work.indices), 4)
    size(data) == size(mask) == shape || throw(DimensionMismatch("expected data shape $shape"))
    sigma = noise_std_mm isa Number ? fill(noise_std_mm, shape) : copy(noise_std_mm)
    size(sigma) == shape || throw(DimensionMismatch("expected noise shape $shape"))
    selected = findall(mask)
    isempty(selected) && throw(ArgumentError("at least one measurement is required"))
    all(i -> isfinite(data[i]) && isfinite(sigma[i]) && sigma[i] > 0, selected) ||
        throw(ArgumentError("invalid observed data or noise"))
    measured = copy(data)
    θ -> begin
        decoded = decode(θ)
        varies_machine = hasproperty(decoded, :machine)
        varies_machine && !isnothing(prepared.fixed_mean) &&
            throw(ArgumentError("cannot vary machine settings with cached transport"))
        predictions = map(experiments) do experiment
            if varies_machine
                set_operating_point!(experiment.work, experiment.op; decoded.machine...)
            end
            experiment.predict(decoded.beam)
        end
        sum(selected) do i
            e, m, c = Tuple(i)
            s = sigma[i]
            residual = (predictions[e][m, c] - measured[i]) / s
            -(residual^2 + log(2π)) / 2 - log(s)
        end
    end
end
