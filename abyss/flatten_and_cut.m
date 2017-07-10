function wave_points_xy = ...
    flatten_and_cut(central_vertex_id, dist, second_on_path, ...
        first_faces, Thetas, bst_surf, max_distance)

    function point = intersection_point(A, B, C, theta)
    % Returns 3D coordinates of point T on side AB such that angle TCA = theta
        BAC = atan2(norm(cross(C-A, B-A)), dot(C-A, B-A));
        CTA = pi - theta - BAC;
        AT = norm(C-A)*sin(theta)/sin(CTA);
        point = A + AT.*(B-A)/norm(B-A);
    end

    %     % For each vertex for which there is a shortest path find the
    % % intersection of that path with the opposite side of the first triangle on
    % % its way. If the path begins with an edge than use the other end of that
    % % edge. We need this intersection points to assing an angle to each point.
    vertex_cnt = size(dist, 1);
    intersections = cell(vertex_cnt, 1);
    for v_id = 1:vertex_cnt
        if dist(v_id) < max_distance
            second_vertex_id = second_on_path(v_id);
            if ~first_faces{central_vertex_id}.isKey(second_vertex_id)
                intersections{v_id} = bst_surf.Vertices(second_vertex_id, :);
                continue;
            end
            first_face = first_faces{central_vertex_id}(second_vertex_id);
            first_face_coord = num2cell(bst_surf.Vertices(first_face, :), 2);
            [A, B, C] = first_face_coord{:};
            theta = Thetas{central_vertex_id}(second_vertex_id);
            intersections{v_id} = intersection_point(A, B, C, theta);
        end
    end
    
%     [central_faces, ~] = find(bst_surf.Faces == central_vertex_id);
%     trimesh(triangulation(bst_surf.Faces(central_faces, :), bst_surf.Vertices));
%     hold on;
%     intersections_locations_mat = cell2mat(intersections);
%     scatter3(...
%         intersections_locations_mat(:, 1), ...
%         intersections_locations_mat(:, 2), ...
%         intersections_locations_mat(:, 3), ...
%         100, 'filled');
%     hold off;
    

    n = vertexNormal( ...
        triangulation(bst_surf.Faces, bst_surf.Vertices), ...
        central_vertex_id);

    flattened_3d_points = containers.Map( ...
            'KeyType', 'int32', 'ValueType', 'any');
    O = bst_surf.Vertices(central_vertex_id, :);
    for v_id = 1:vertex_cnt
        if dist(v_id) < max_distance && v_id ~= central_vertex_id
            T = intersections{v_id};
            P = T - dot(n, T-O)*n;
            flattened_3d_points(v_id) = O + (P-O)/norm(P-O)*dist(v_id);
        end
    end
    flattened_3d_points(central_vertex_id) = O;

% animation of flattening
%     for i = 0:10
%         vertices = bst_surf.Vertices;
%         vertices([wave_vertices_idx{:}], :) = ...
%             vertices([wave_vertices_idx{:}], :) .* ((10-i)/10) ...
%             + cell2mat(flattened_3d_points').* (i/10);
%         faces = bst_surf.Faces;
%         
%         index_mapping_anim = 1:vertex_cnt;
%         index_mapping_anim([wave_vertices_idx{:}]) ...
%             = 1:length(wave_vertices_idx);
%         index_mapping_anim(setdiff(1:vertex_cnt, [wave_vertices_idx{:}])) = NaN; 
% 
%         vertices = vertices([wave_vertices_idx{:}], :);
%         faces = arrayfun(@(old_vertex_id) index_mapping_anim(old_vertex_id), ...
%             faces);
%         faces(any(isnan(faces), 2), :) = [];
%         
%         
%         trimesh(triangulation(faces, vertices)); 
%         pause(eps);
%     end

    % Move to 2D
    wave_vertices_idx = cell2mat(flattened_3d_points.keys);
    random_wave_vertex_idx = wave_vertices_idx( ...
        find(wave_vertices_idx ~= central_vertex_id, 1));
    random_flattend_point = flattened_3d_points(random_wave_vertex_idx);
    x = random_flattend_point - O;
    x = x/norm(x);
    y = cross(n, x);
    
    wave_points_xy = containers.Map( ...
            'KeyType', 'int32', 'ValueType', 'any');
    for key = flattened_3d_points.keys()
        wave_points_xy(key{1}) = [x; y] * (flattened_3d_points(key{1}) - O)';
    end

end