function betas = ...
    lasso_the_waves(sensor_waves, wave_inds, ...
    start_times, meas, PARAMS)

    equalize_timewise = PARAMS.lasso.equalize_timewise;
    time_point_cnt = get_time_point_cnt(PARAMS);
    zero_point_id = find(abs(meas{1}.Time) < 10e-6);
    
    sensor_waves_sz = size(sensor_waves);
    sensor_waves = reshape(sensor_waves, prod(sensor_waves_sz(1:3)), ...
        sensor_waves_sz(4), sensor_waves_sz(5));
    sensor_waves = sensor_waves(wave_inds, :, :);
    sensor_waves = padarray(sensor_waves, [0 0 time_point_cnt-1], 0, 'post');
    
    trial_cnt = size(meas, 2);
    betas = cell(trial_cnt, 1);
    wave_cnt = size(sensor_waves, 1);
    for t_id = 1:trial_cnt
        
        fprintf('Final WLS for trial %d/%d\n', t_id, trial_cnt);
        
        y = meas{t_id}.F(:, ...
            zero_point_id + [-time_point_cnt+1 : time_point_cnt-1]);
        X = sensor_waves;
        st = start_times(:, t_id);
        for w_id = 1:wave_cnt
            X(w_id, :, :) = circshift(X(w_id, :, :), st(w_id)-1, 3);
        end
        if equalize_timewise
            sigma = sqrt(sum(y.^2, 1));
            y = bsxfun(@rdivide, y, sigma);
            X = bsxfun(@rdivide, X, shiftdim(sigma, -1));
        end
        y = y(:);
        X = reshape(X, wave_cnt, [])';
        
        % OLS
        betas{t_id} = X\y;
        
    end
end