function states = KalmanFilter(t, signal, noise_std, q)
    n = length(t);
    
    dim = size(signal, 1);

    err_measure = noise_std;
    err_estimate = err_measure;

    states = zeros(size(signal));

    last_estimate = signal(:, 1);

    for i = 1:n
        mea = signal(:, i);
        
        E = err_measure + err_estimate;

        kalman_gain = err_estimate ./ E;
        K = diag(kalman_gain);

        curr_estimate = last_estimate + K * (mea - last_estimate);

        diff = diag(abs(last_estimate - curr_estimate));
        err_estimate = (diag(ones(dim,1)) - K)*err_estimate + diff*q;

        last_estimate = curr_estimate;
        states(:, i) = curr_estimate;
    end
end