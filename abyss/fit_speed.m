function speed_fitted = fit_speed(planar_mesh, meas, gain_matrix, PARAMS)
% Find a speed of propagation for which the wave going in all directions
% fits the measurements meas the best

    half_width = PARAMS.wave_defaults.half_width;
    direction = NaN; angular_half_width = NaN;
    % Number of column minus 1 is the number of sampling intervals. Minus 1
    % because the first column corresponds to time 0.
    time_interval_ms = (size(meas, 2)-1)/PARAMS.sampling_rate*1000;
    generate_wave_with_given_speeed = @(speed) ...
        generate_basis_wave( ...
                planar_mesh.Vertices, planar_mesh.Faces, ...
                time_interval_ms, speed, half_width, ...
                direction, angular_half_width, ...
                PARAMS.sampling_rate);
            
    speeds = 0.1:0.1:2;
    waves_with_various_speeds = {};
    for speed = speeds
        try
            waves_with_various_speeds{end+1} = ...
                generate_wave_with_given_speeed(speed);
        catch ME
            if (strcmp(ME.identifier, 'generate_basis_wave:reachedBorder'))
                break;
            else
                rethrow(ME);
            end
        end
    end
        

    speed_coeff = pick_waves_using_lasso(waves_with_various_speeds, meas, ...
        gain_matrix, []);

    [~, speed_fitted_id] = max(speed_coeff); 
    speed_fitted = speeds(speed_fitted_id);
    
    if ismember(speed_fitted_id, [1, length(waves_with_various_speeds)])
        warning(['Speed fitted (%0.2f mm/ms) is either the slowest or the ' ...
            'fastest'], speed_fitted);
    end
end