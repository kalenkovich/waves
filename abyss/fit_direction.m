function [wave_fitted, fit_info] = fit_direction(planar_mesh, meas, ...
    gain_matrix, speed, PARAMS)

    half_width = PARAMS.wave_defaults.half_width;
    angular_half_width = PARAMS.wave_defaults.angular_half_width;
    % Number of column minus 1 is the number of sampling intervals. Minus 1
    % because the first column corresponds to time 0.
    time_interval_ms = (size(meas, 2)-1)/PARAMS.sampling_rate*1000;
    generate_wave_with_given_direction = @(direction) ...
        generate_basis_wave( ...
                planar_mesh.Vertices, planar_mesh.Faces, ...
                time_interval_ms, speed, half_width, ...
                direction, angular_half_width, ...
                PARAMS.sampling_rate);

    % generate waves in different directions
    direction_cnt = round(2*pi/angular_half_width);
    directions = 2*pi/direction_cnt .* [0:direction_cnt-1];
    directions(end+1) = NaN; % NaN means the propagation is in all directions

    waves = arrayfun(generate_wave_with_given_direction, directions, 'un', 0);

    [coeff, fit_info] = pick_waves_using_lasso(waves, meas, gain_matrix, ...
        direction_cnt + 1);
    
    % Check that only one direction is picked
    direction_id = find(coeff(1:end-1));
    if length(direction_id) ~= 1
        warning('More than one direction picked. You''d better investigate');
    end

    % Multiply each wave by its coefficient, then add them all up
    wave_fitted = arrayfun(@(wave, q) wave{:}*q, waves, coeff', 'un', 0);
    wave_fitted = sum(cat(3, wave_fitted{:}), 3);

    % Find the vertex reached by the directed wave with the highest
    % coefficient. To be exact - the vertex closest to the peak of that
    % wave.
    [~, direction_id] = max(coeff(1:end-1));
    direction = directions(direction_id);
    last_time_point = roundn((size(meas, 2)-1) / PARAMS.sampling_rate, -3);
    actual_peak = speed * last_time_point ...
        .* [cos(direction), sin(direction)];
    [~, reached_vertex_id] = min( ...
        sum(bsxfun(@minus, planar_mesh.Vertices, actual_peak).^2, 2));
    
    % Collect information into fit_info
    fit_info.reached_vertex_id = reached_vertex_id;
    
end