function wave = generate_basis_wave(vertices_xy, faces, time_interval_in_ms, ...
    speed, half_width, direction, angular_half_width, SAMPLING_RATE)
%     Returns a matrix in which each row corresponds to one vertex and each
%     column - to a time point. 
%
%     @param wave_points_xy : map of vertex indices and their planar coordinates
%     @param wave_faces : triangles of the planar mesh given as indices of each vertex in the points array
%     time : an int array of T time points in time in seconds
%     speed : an int speed of the wave crest in mm/ms
%     half_width : the wave is non-zero for points such that ||point| - time*speed| < half_width
%     
%     direction : optional angle in radians corresponding to the middle of the crest's arc
%     angular_half_width : optional angular half-width of the crest's arc. 
%    
%     Wave as a raised cosine in each radial direction centered on the 
%     crest and supported in its neigbourhood with halfwidth as radius
    
    point_cnt = size(vertices_xy, 1);
    sample_cnt = ceil(SAMPLING_RATE/1000 * time_interval_in_ms) + 1;
    time_points = 1/SAMPLING_RATE * [0:sample_cnt-1];
    
    % @var mu and r are (point_cnt X sample_cnt)-size array
    % @var mu contains distances from zero to the crest for each time point
    % @var r contains distance from zero to each vertex
    
    
    mu = time_points*speed;
    mu = repmat(mu, point_cnt, 1);
    
    wave_point_idx = 1:point_cnt;
    r = sqrt(sum(vertices_xy.^2, 2));
    r = repmat(r, 1, sample_cnt);

    raised_cosine = @(x, mu, s) (1 + cos((x-mu) / s*pi)) ...
        .* (mu-s < x) .* (x < mu+s);
    
    wave = ... 
        raised_cosine(r, mu, half_width) ... 
        + raised_cosine(r, -mu, half_width);
    

%      If it is a directed wave we apply a mask that is agained a raised 
%      cosine but now on any circle centered at zero, the center of the
%      raised cosine corresponds to the point on the circle whose angle is
%      direction of the wave and the support is the arc that is intercepted
%      by angular neigbourhood of direction with radius angular_half_width.
%      Center point is a special case and should not be masked.

    angle_of = @(point) atan2(point(2), point(1));    
    if ~isnan(direction) && ~isnan(angular_half_width)
        phi = arrayfun(@(key) angle_of(vertices_xy(key, :)), ...
            wave_point_idx);
        phi(find(r(:,1) == 0)) = direction;
        angle_deltas = mod((phi-direction+pi), 2*pi) - pi; % in [-pi,pi)
        direction_mask = raised_cosine(angle_deltas, 0, angular_half_width);
        wave = wave.* repmat(direction_mask', 1, sample_cnt);
    end

%      Now we normalize the waves so that instanstaneous energy (defined here
%      as the integral of the square of the wave) is equal to zero. We
%      calculate the integral as the volume under the wave mesh
    
    % Coordinate of all the vertices of the faces of the wave
    % tri_xy: {N triangles x 3 vertices} x [2 coordinates]
    tri_xy = arrayfun(@(p_id) vertices_xy(p_id, :), faces, 'un', 0);
    % tri_z: {N triangles x 3 vertices} x [1 coordinate x T time points]
    tri_z = arrayfun(@(p_id) wave(p_id, :), faces, 'un', 0);
    % tri_z: [N triangles x 3 vertices x T time points]
    tri_z = reshape(cat(1, tri_z{:}), [size(tri_z), sample_cnt]);

    % For each face there is a prism. We need to calculate its average
    % heigth and the area of its base.
    average_height = squeeze(mean(tri_z.^2, 2)); % [N triangles x T time points]
    a = arrayfun(@(i) tri_xy{i, 2} - tri_xy{i, 1}, 1:size(tri_xy, 1), 'un', 0);
    b = arrayfun(@(i) tri_xy{i, 3} - tri_xy{i, 1}, 1:size(tri_xy, 1), 'un', 0);
    tri_area = @(a, b) 0.5 * abs(a(1)*b(2) - a(2)*b(1));
    base_area = cellfun(tri_area, a, b)'; % [N triangles x 1] 
    
    % Summing the volumes of all the prisms we find the integral
    energy = sum(bsxfun(@times, average_height, base_area)); % [1 x T time points]
    wave = bsxfun(@rdivide, wave, sqrt(energy)); % [M vertices x T time points]
    
    % Check that the mesh is wide enough, i.e. that the wave does not reach
    % the border
    TR = triangulation(faces, vertices_xy);
    boundary_vertices_idx = unique(TR.freeBoundary());
    boundary_values = wave(boundary_vertices_idx, :);
    if any(boundary_values(:) ~= 0)
        errorStruct.message = ['The wave generated reached the boundary ' ...
            'which is bad in current implementation because of the way the waves ' ...
            'are normalized. Either widen the submesh on which the waves are generated ' ...
            'or lower the speed.'];
        errorStruct.identifier = 'generate_basis_wave:reachedBorder';
        error(errorStruct);
    end
        

 end