function avg_orientation = averageOrientationsLogEuclidean(orientations)

    eps = 1e-6;
    log_avg_tensor = zeros(3);
    nb_orientations = size(orientations, 1);
    for k = 1:nb_orientations
        current_tensor = orientations(k, :)' * orientations(k, :) + eps*eye(3);
        log_avg_tensor = log_avg_tensor + logm(current_tensor);
    end
    avg_tensor = expm(log_avg_tensor / nb_orientations);
    
    [V, D] = eigs(avg_tensor);
    avg_orientation = V(:, 1)';
end