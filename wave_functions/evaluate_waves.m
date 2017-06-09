function scores = evaluate_waves(sensor_waves, meas, PARAMS)

    % Extract params
    speeds = PARAMS.wave.speeds;
    speed_cnt = length(speeds);
    directions = PARAMS.wave.directions;
    direction_cnt = length(directions);
    time_point_cnt = get_time_point_cnt(PARAMS);
    meas_cnt = length(meas);
    
    intercept = PARAMS.scores.intercept;
    equalize_timewise = PARAMS.scores.equalize_timewise;
    
    scores = cell(meas_cnt, 1);
    % Loop over measurements in meas
    for m_id = 1:meas_cnt
    
        % Cut the measurements with a sliding window of the same width as the
        % duration of waves. Do not start later than 0 though.
        zero_point_id = find(abs(meas{m_id}.Time) < 10e-6);
        meas_windowed = arrayfun(@(t) meas{m_id}.F(:, t:t+time_point_cnt-1), ...
                         1:zero_point_id, 'un', 0);
        meas_windowed = squeeze(cell2ndim(meas_windowed));

        % Remove mean from each window and each wave
        if intercept
            meas_means = mean(mean(meas_windowed, 2), 3);
            meas_windowed = bsxfun(@minus, meas_windowed, meas_means);
            wave_means = mean(mean(sensor_waves, 4), 5);
            sensor_waves = bsxfun(@minus, sensor_waves, wave_means);
        end

        % Divide at each time point by the norm of the measurement at this time
        % point. We cannot apply weights to the waves at this point because the
        % weights depend on the window.
        window_cnt = size(meas_windowed, 1);
        time_dim = ndims(meas_windowed);
        sensor_dim = time_dim - 1;
        if equalize_timewise
            weight_by_time = sqrt(sum(meas_windowed.^2, sensor_dim)); 
        else
            weight_by_time = ones(window_cnt, 1, time_point_cnt);
        end
        meas_windowed = bsxfun(@rdivide, meas_windowed, weight_by_time);

        % Convert to a matrix where each window is one row
        meas_windowed = reshape(meas_windowed, window_cnt, []);

        % Normalize the measurements so that the norm in each window is 1.
        meas_norms = sqrt(sum(meas_windowed.^2, 2));
        meas_windowed = bsxfun(@rdivide, meas_windowed, meas_norms);

        % Calculate the scores. There is one score for each wave, time window
        % [and measurement]

        sensor_waves_sz = size(sensor_waves);
        scores{m_id} = nan([sensor_waves_sz(1:end-2) window_cnt]);

        time_dim_of_waves = ndims(sensor_waves);
        sensor_cnt = sensor_waves_sz(end-1);
        for w_id = 1:window_cnt
            tic;
            fprintf('Evaluating waves against the measurements %d/%d on window %d/%d\n', ...
                m_id, meas_cnt, w_id, window_cnt);

            % Apply weights
            w = weight_by_time(w_id, :)';
            w = shiftdim(w, 1 - time_dim_of_waves);
            waves = bsxfun(@rdivide, sensor_waves, w);

            % Convert to a matrix where each waves is one row
            waves_sz = [prod(sensor_waves_sz(1:end-2)), ...
                        sensor_cnt * time_point_cnt];
            waves = reshape(waves, waves_sz);

            % Normalize
            wave_norms = sqrt(sum(waves.^2, 2));
            waves = bsxfun(@rdivide, waves, wave_norms);

            % The scores are now the square of the dot product of waves and the
            % measurements in current window
            scores{m_id}(:, :, :, w_id) = reshape(waves * meas_windowed(w_id, :)', ...
                sensor_waves_sz(1:3));
            toc;
        end
    end
end