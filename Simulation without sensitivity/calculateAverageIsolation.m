function avgIsolation = calculateAverageIsolation(iso, nonIso)
    % calculateAverageIsolation computes the average isolation ratio across all frequencies.
    %
    % Inputs:
    %   iso    - Vector of isolated movement data
    %   nonIso - Vector of non-isolated movement data
    %   fs     - Sampling frequency (Hz)
    %
    % Output:
    %   avgIsolation - Weighted average isolation ratio (0 = full isolation, 1 = no isolation)

    % Input validation
    if length(iso) ~= length(nonIso)
        error('The iso and nonIso vectors must be of the same length.');
    end

    N = length(iso);  % Signal length

    % Handle DC component (mean offset)
    if all(iso == iso(1))
        % iso is constant
        avgIsolation = mean(iso) / mean(nonIso); % DC ratio only
    else
        % Remove DC offset for dynamic analysis
        iso = iso - mean(iso);
        nonIso = nonIso - mean(nonIso);

        % FFT of signals
        fftIso = abs(fft(iso)) / N;
        fftNonIso = abs(fft(nonIso)) / N;

        % Handle division by zero
        zeroIdx = fftNonIso == 0;
        fftNonIso(zeroIdx) = NaN;  % Prevent division by zero

        % Compute isolation ratio
        isolationRatio = fftIso ./ fftNonIso;

        % Restore DC ratio (mean comparison at 0 Hz)
        isolationRatio(1) = mean(iso + mean(iso)) / mean(nonIso + mean(nonIso));

        % Compute weighted average isolation ratio
        validIdx = ~isnan(isolationRatio);
        avgIsolation = sum(isolationRatio(validIdx) .* fftNonIso(validIdx)) / sum(fftNonIso(validIdx));
    end
end
