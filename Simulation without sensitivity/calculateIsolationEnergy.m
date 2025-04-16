function isolationPercent = calculateIsolationEnergy(iso, nonIso)
    % Calculates the percent of original signal energy retained in the controlled signal.
    % 0% = full isolation (iso is 0), 100% = no isolation (iso == nonIso)

    if length(iso) ~= length(nonIso)
        error('Vectors must be the same length.');
    end

    % Remove mean (DC offset)
    iso = iso - mean(iso);
    nonIso = nonIso - mean(nonIso);

    % Compute signal energies
    energyIso = sum(iso.^2);
    energyNonIso = sum(nonIso.^2);

    % Ratio of energies
    isolationPercent = energyIso / energyNonIso;  % 1 = no isolation
end
