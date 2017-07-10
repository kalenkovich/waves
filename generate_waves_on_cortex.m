function cortex_waves = generate_waves_on_cortex(cortex, PARAMS)
% Generates waves on the cortex and project them onto sensors. At each time
% point the waves are normalized by the square root of the squared wave
% integrated over the surface. This integral is calculated on a wave
% generated on a plane.


    % Load the forward-model matrix
    G = from_bst_get_gain_matrix(PARAMS.forward.name, PARAMS);

    % Zeros on the diagonal of the distance matrix represent actual zeros
    % while the other zeros are actually infinities. For angles the scheme
    % has additional twists. Thus it is easier to work with the (row index,
    % column index, value) triplets.
    [i, j, r] = find(cortex.distances);
    phi = full(cortex.angles(cortex.distances ~= 0));
    vertex_cnt = size(cortex.distances, 1);
    i = vertcat(i, [1:vertex_cnt]');
    j = vertcat(j, [1:vertex_cnt]');
    r = vertcat(r, zeros(vertex_cnt, 1));
    phi = vertcat(phi, zeros(vertex_cnt, 1));
    
    % This is the meshgrid we will generate waves on for normalization
    % purposes
    max_distance = PARAMS.max_distance;
    step = 0.0001;
    x = -max_distance:step:max_distance;
    y = -max_distance:step:max_distance;
    [X, Y] = meshgrid(x, y);
    assert(nnz(X==0 & Y==0) == 1, ['Meshgrid used to calculate norms of ' ...
        'the waves does not include the center. Probably weird max_distance.']);
    [norm_phi, norm_r] = cart2pol(X, Y);
    
    % Extract parameters from the PARAMS struct
    duration = PARAMS.wave.duration;
    half_width = PARAMS.wave.half_width;
    angular_half_width = PARAMS.wave.angular_half_width;
    speeds = PARAMS.wave.speeds;
    directions = PARAMS.wave.directions;
    
    time_points = 0:1/PARAMS.sampling_rate:duration;
    
   
    
    % When we calculate the direction mask the 0 angle that corresponds to
    % the center point is a special case and we assing the direction angle
    % to it. For the cortex these corresponds to points where i==j, for the
    % mesh grid that we use to normalize the wave this is the point where
    % Y==0 & X==0 or equivalently (Y==0 & X==0) == 1.
    direction_masks = arrayfun(@(dir) ...
        generate_direction_mask(i, j, phi, dir, angular_half_width), ...
        directions, 'un', 0);
    norm_direction_masks = arrayfun(@(dir) ...
        generate_direction_mask(1, Y==0 & X==0, ...
                                norm_phi, dir, angular_half_width), ...
        directions, 'un', 0);
    
    % Empty cell array to hold the results
    time_point_cnt = length(time_points);
    speed_cnt = length(speeds);
    direction_cnt = length(directions);
    cortex_waves = cell(speed_cnt, direction_cnt, time_point_cnt);
    
    % Main loop
    for speed_id = 1:speed_cnt
        for time_point_id = 1:time_point_cnt
            disp(sprintf('Generating waves for speed %d and time %d', ...
                speed_id, time_point_id));
            
            speed = speeds(speed_id);
            t = time_points(time_point_id);
            
            ripple_wave_v = generate_ripple_wave(r, speed, t, half_width);
            norm_ripple_wave = ...
                generate_ripple_wave(norm_r, speed, t, half_width);
            
            norm_waves = cellfun(@(norm_dir_mask) ...
                norm_ripple_wave .* norm_dir_mask, norm_direction_masks, ...
                'un', 0);
            integral_over_grid = @(x, y, z) trapz(y, trapz(x, z, 2));
            norm_factors = cellfun(@(norm_wave) ...
                sqrt(integral_over_grid(x, y, norm_wave)), norm_waves, 'un', 0);
            
            waves = cellfun(@(dir_mask, factor) ...
                sparse(i, j, ripple_wave_v .* dir_mask / factor), ...
                direction_masks, norm_factors, 'un', 0);
            
            [cortex_waves{speed_id, :, time_point_id}] = ...
                deal(waves{:});
            
        end
    end
    
    % Right night it goes {speed, direction, time_point}[vertex, sensor]
    % which does not make a lot of sense. Let's make it right by converting
    % this monstrosity into [vertex, speed, direction, sensor, time_point]
    % multi-dimensional array.
    % TODO: change the function so that sensor_waves is this way from the
    % start.
%     cortex_waves = cellfun(@(x) reshape(x, 1, 1, 1, size(x, 1), size(x, 2)), ...
%         cortex_waves, 'un', 0);
%     cortex_waves = cell2mat(cortex_waves);
%     cortex_waves = permute(cortex_waves, [4, 1, 2, 5, 3]);
    
end