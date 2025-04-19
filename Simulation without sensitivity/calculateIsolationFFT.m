function isolationPercent = calculateIsolationFFT(iso, nonIso)
    if length(iso) ~= length(nonIso)
        error('Vectors must be the same length.');
    end

    N = length(iso);
    iso = iso - mean(iso);
    nonIso = nonIso - mean(nonIso);

    % FFT
    fftIso = abs(fft(iso)/N).^2;
    fftNonIso = abs(fft(nonIso)/N).^2;

    % Single-sided spectrum
    fftIso = fftIso(1:floor(N/2)+1);
    fftNonIso = fftNonIso(1:floor(N/2)+1);

    % Energy ratio
    isolationPercent = sqrt(sum(fftIso) / sum(fftNonIso));  % 1 = no isolation
end
