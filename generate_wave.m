function wave = generate_wave(central_vertex_id, speed_id, direction_id, ...
    cortex, PARAMS)

    % Extract parameters from the PARAMS struct
    v_id = central_vertex_id;
    speed = PARAMS.wave.speeds(speed_id);
    direction = PARAMS.wave.directions(direction_id);
    duration = PARAMS.wave.duration;
    half_width = PARAMS.wave.half_width;
    angular_half_width = PARAMS.wave.angular_half_width;
    
    time_points = 0:1/PARAMS.sampling_rate:duration;
    
    % We only need distances and angles from one vertex
    cortex.distances = cortex.distances(v_id, :);
    cortex.angles = cortex.angles(v_id, :);
    
    % Zeros on the diagonal of the distance matrix represent actual zeros
    % while the other zeros are actually infinities. For angles the scheme
    % has additional twists. Thus it is easier to work with the (row index,
    % column index, value) triplets.
    [i, j, r] = find(cortex.distances);
    phi = full(cortex.angles(cortex.distances ~= 0));
    vertex_cnt = size(cortex.distances, 2);
    % There is only one row in cortex.distance and cortex.angles therefore
    % i, j, r and phi are row vectors. Fuck MATLAB!
    i = [i'; 1];
    j = [j'; v_id];
    r = [r'; 0];
    phi = [phi'; 0];
    
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
    
    % When we calculate the direction mask the 0 angle that corresponds to
    % the center point is a special case and we assing the direction angle
    % to it. For the cortex these corresponds to points where i==j, for the
    % mesh grid that we use to normalize the wave this is the point where
    % Y==0 & X==0 or equivalently (Y==0 & X==0) == 1.
    %
    % The above is true when we generate waves all over the cortex. Here
    % the central point corresponds to j == v_id.
    direction_mask = generate_direction_mask(v_id, j, phi, ...
        direction, angular_half_width);
    norm_direction_mask = generate_direction_mask(1, Y==0 & X==0, norm_phi, ...
        direction, angular_half_width);
    
    time_point_cnt = get_time_point_cnt(PARAMS);
    wave = sparse(vertex_cnt, time_point_cnt);
    for time_point_id = 1:time_point_cnt

        t = time_points(time_point_id);

        ripple_wave_v = generate_ripple_wave(r, speed, t, half_width);
        norm_ripple_wave = ...
            generate_ripple_wave(norm_r, speed, t, half_width);

        norm_wave = norm_ripple_wave .* norm_direction_mask;
        integral_over_grid = @(x, y, z) trapz(y, trapz(x, z, 2));
        norm_factor = sqrt(integral_over_grid(x, y, norm_wave));

        wave(:, time_point_id) = sparse(i, j, ...
            ripple_wave_v .* direction_mask / norm_factor, ...
            1, vertex_cnt);

    end
    
    
    
end