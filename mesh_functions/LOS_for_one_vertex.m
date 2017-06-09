function [min_geod_distances, first_triangles, first_angles] = ...
    LOS_for_one_vertex(mesh, S_id, max_distance)
% Calculates minimal-geodesic distances from the vertex with index v_id not 
% longer than max_distance.
% Also retains information about the first triangle that the path goes
% through in the form (for one destination vertex) [S_id, A_id, B_id] and
% angles between SA and the first part of the path.

    % Plot for debugging
    % DEBUG: fig = figure; hold on;

    % Array allocation
    vertex_cnt = length(mesh.Vertices);
    min_geod_distances = inf(vertex_cnt, 1);
    first_triangles = zeros(vertex_cnt, 3);
    first_angles = NaN(vertex_cnt, 1);

    % Initialze some of the variables
    vertex_template = struct('id', NaN, 'xyz', NaN, 'xy', NaN);
    [A, B, C, D, S] = deal(vertex_template);
    S.id = S_id; S.xyz = mesh.Vertices(S_id, :); S.xy = [0, 0];
    min_geod_distances(S.id) = 0;
    stack = {};
    stack_size = 0;
    TR = triangulation(mesh.Faces, mesh.Vertices);
    
    % Precalculate some values
    mesh.angles = precalculate_angles();
    mesh.Faces_sorted = sort(mesh.Faces, 2);
    [first_attached_vertex, second_attached_vertex] = ...
        precalculate_attached_vertices();
    
    % Main loop
    adjacent_faces_idx = get_adjacent_faces_idx(S_id);    
    adjacent_faces_cnt = length(adjacent_faces_idx);
    for Tau_id = 1:adjacent_faces_cnt
        face_id = adjacent_faces_idx(Tau_id);
        Tau = mesh.Faces(face_id, :);
        [first_triangle, theta_max] = process_first_triangle();
        theta_min = 0;
        push_to_stack(S, A, B, theta_min, theta_max, face_id);
        while stack_size > 0
            pop_from_stack();
            % DEBUG: disp(num2str([C.id A.id B.id]));
            [D, theta] = extend_over_AB();
            if isnan(D.id) % || (D.id == S.id)
                continue;
            end
            if theta <= theta_max
                push_to_stack(A, D, B, max(theta, theta_min), theta_max, face_id);
            end
            if theta >= theta_min
                push_to_stack(B, A, D, theta_min, min(theta, theta_max), face_id);
            end
            if theta <= theta_max && theta >= theta_min
                update_info_about_D();
            end
        end
    end
    
    % There is no first_triangle or angle for S.id, so here I remove min
    % distance for S.id for consistency of outputs.
    min_geod_distances(S.id) = inf;

    function [first_triangle, theta_max] = process_first_triangle()
        
        % Fill A and B
        AB = Tau(Tau ~= S_id);
        A.id = AB(1); 
        B.id = AB(2);
        A.xyz = mesh.Vertices(A.id, :);
        B.xyz = mesh.Vertices(B.id, :);
        A.xy = [norm(A.xyz - S.xyz), 0];
        B.xy = map_triangle_to_the_plane(A, S, B);
        
        % Save Tau in a format that we need for first_triangles
        first_triangle = [S.id, A.id, B.id];
        theta_max = atan2(B.xy(2), B.xy(1));
        
        % Update info on paths to A and B
        min_geod_distances(A.id) = norm(A.xyz - S.xyz);
        min_geod_distances(B.id) = norm(B.xyz - S.xyz);
        first_triangles(A.id, :) = first_triangle;
        first_triangles(B.id, :) = first_triangle;
        first_angles(A.id) = 0;
        first_angles(B.id) = theta_max;
    end

    function C_xy = map_triangle_to_the_plane(A, B, C)
    % Returns planar coordianates for the point C such that: 1) triangles
    % with planar and with spatial coordinates of A, B and C are congruent;
    % 2) C.xy is on the right of A.xy -> B.xy

        % ABC = angle_between_3d_vectors(A.xyz - B.xyz, C.xyz - B.xyz);
        % ABC = get_angle_in_face(face_id, B.id);
        ABC = mesh.angles((A.id-1)*vertex_cnt + B.id, C.id);
        c = cos(ABC); s = sin(ABC);
        R = [c, -s; s, c]; % Rotation matrix
        BA_xy_unit = (A.xy-B.xy)' ./ norm(A.xy-B.xy); % Planar unit vector in direction BA
        C_xy = B.xy' + R * BA_xy_unit .* norm(C.xyz-B.xyz);
        C_xy = C_xy';
    end

    function angle = angle_between_3d_vectors(a, b)
        % Copy-pasted from https://goo.gl/d9DTcC
        % angle = atan2(norm(cross(a, b)), dot(a, b));
        angle = 2 * atan( ...
            norm(a*norm(b) - norm(a)*b) ...
          / norm(a * norm(b) + norm(a) * b));
    end

    function angles = precalculate_angles()
        % for each vertex of each face calculate a and b - side of the
        % faces as vectors from this vertex
        a = cat(3, ...
            mesh.Vertices(mesh.Faces(:, 2), :) ...
                - mesh.Vertices(mesh.Faces(:, 1), :), ...
            mesh.Vertices(mesh.Faces(:, 3), :) ...
                - mesh.Vertices(mesh.Faces(:, 2), :), ...
            mesh.Vertices(mesh.Faces(:, 1), :) ...
                - mesh.Vertices(mesh.Faces(:, 3), :));
        b = cat(3, ...
            mesh.Vertices(mesh.Faces(:, 3), :) ...
                - mesh.Vertices(mesh.Faces(:, 1), :), ...
            mesh.Vertices(mesh.Faces(:, 1), :) ...
                - mesh.Vertices(mesh.Faces(:, 2), :), ...
            mesh.Vertices(mesh.Faces(:, 2), :) ...
                - mesh.Vertices(mesh.Faces(:, 3), :));
            
        % Apply permutations so that the dimensions are face_id, vertex, axis
        a = permute(a, [1 3 2]);
        b = permute(b, [1 3 2]);

        % Apply formula copy-pasted from https://goo.gl/d9DTcC
        a_norm = sqrt(sum(a.^2, 3));
        b_norm = sqrt(sum(b.^2, 3));

        axb_norm = bsxfun(@times, a, b_norm);
        bxa_norm = bsxfun(@times, b, a_norm);
        
        axb_norm_minus_bxa_norm = axb_norm - bxa_norm;
        axb_norm_plus_bxa_norm = axb_norm + bxa_norm;
        
        minus_norm = sqrt(sum(axb_norm_minus_bxa_norm.^2, 3));
        plus_norm = sqrt(sum(axb_norm_plus_bxa_norm.^2, 3));
        
        ratio = minus_norm ./ plus_norm;
        angles = 2 * atan(ratio); 
        
        % Create a sparse matrix that has angle (A_id, B_id, C_id) at
        % position ((A_id - 1) * vertex_cnt + B_id, C_id) - kind of a 3D
        % sparse matrix
        A_ids = vertcat(...
            mesh.Faces(:, 3), mesh.Faces(:, 1), mesh.Faces(:, 2), ...
            mesh.Faces(:, 2), mesh.Faces(:, 3), mesh.Faces(:, 1));
        B_ids = vertcat(...
            mesh.Faces(:, 1), mesh.Faces(:, 2), mesh.Faces(:, 3), ...
            mesh.Faces(:, 1), mesh.Faces(:, 2), mesh.Faces(:, 3));
        C_ids = vertcat(...
            mesh.Faces(:, 2), mesh.Faces(:, 3), mesh.Faces(:, 1), ...
            mesh.Faces(:, 3), mesh.Faces(:, 1), mesh.Faces(:, 2));
        angles = vertcat(...
            angles(:, 1), angles(:, 2), angles(:, 3), ...
            angles(:, 1), angles(:, 2), angles(:, 3));
        AB_ids = (A_ids - 1) * vertex_cnt + B_ids;
        angles = sparse(AB_ids, C_ids, angles);
    end

    function [first_attached_vertex, second_attached_vertex] = ...
            precalculate_attached_vertices()
    % For each edge finds its adjacent faces and saves third vertices of
    % those faces in two sparse matrices. Second one has less elements.
        
        mesh_edges = edges(TR);
        attachments = edgeAttachments(TR, mesh_edges);
        first_faces = cell2mat(cellfun(...
            @(att) mesh.Faces(att(1), :), attachments, 'un', 0));
        
        attachment_cnts = cellfun(@length, attachments);
        assert(max(attachment_cnts) == 2, ...
            'Handles in the mesh. Damn you brainstorm!');
        
        edges_with_two_faces = mesh_edges(attachment_cnts == 2, :);
        two_attachments = attachments(attachment_cnts == 2);
        second_faces = cell2mat(cellfun(...
            @(att) mesh.Faces(att(2), :), two_attachments, 'un', 0));
        
        % From each row remove vertices of the corresponing edge
        first_vertex = sum(first_faces .* ...
             bsxfun(@ne, first_faces, mesh_edges(:, 1)) .*...
             bsxfun(@ne, first_faces, mesh_edges(:, 2)), 2);
        second_vertex = sum(second_faces .* ...
             bsxfun(@ne, second_faces, edges_with_two_faces(:, 1)) .*...
             bsxfun(@ne, second_faces, edges_with_two_faces(:, 2)), 2);    
        
        % Now, make sparse matrices from this and account for two possible
        % orderings
        first_attached_vertex = sparse(...
            vertcat(mesh_edges(:, 1), mesh_edges(:, 2)), ...
            vertcat(mesh_edges(:, 2), mesh_edges(:, 1)), ...
            vertcat(first_vertex, first_vertex), ...
            vertex_cnt, vertex_cnt);
        second_attached_vertex = sparse(...
            vertcat(edges_with_two_faces(:, 1), edges_with_two_faces(:, 2)), ...
            vertcat(edges_with_two_faces(:, 2), edges_with_two_faces(:, 1)), ...
            vertcat(second_vertex, second_vertex), ...
            vertex_cnt, vertex_cnt);
    end

    function push_to_stack(C, A, B, theta_min, theta_max, face_id)
        if abs(theta_max - theta_min) > eps(theta_min)*4
            stack_size = stack_size + 1;
            stack{stack_size} = {C, A, B, theta_min, theta_max, face_id}; 
            % DEBUG: triplot([1 2 3], [C.xy(1), A.xy(1), B.xy(1)]', [C.xy(2), A.xy(2), B.xy(2)]');
        end
    end

    function pop_from_stack()
        [C, A, B, theta_min, theta_max, face_id] = stack{stack_size}{:};
        stack_size = stack_size - 1;
    end

    function [D, theta] = extend_over_AB()
        if (norm(A.xy) > max_distance) && (norm(B.xy) > max_distance) ...
                && (distance_to_segment(A, B) > max_distance)
            D.id = NaN;
            theta = NaN;
            return;
        end 
        [C_and_D] = nonzeros([first_attached_vertex(A.id, B.id) ...
            second_attached_vertex(A.id, B.id)]);
        if length(C_and_D) == 2
            D.id = C_and_D(C_and_D ~= C.id);
            D.xyz = mesh.Vertices(D.id, :);
            D.xy = map_triangle_to_the_plane(A, B, D);
            theta = atan2(D.xy(2), D.xy(1));
        elseif length(C_and_D) == 1
            D.id = NaN;
            theta = NaN;
        end           
    end

    function distance = distance_to_segment(A, B)
        AB = B.xy - A.xy;
        AB_norm = norm(AB);
        AB_unit = AB./AB_norm;
        AP_coord = dot(AB_unit, -A.xy); % coordinate on (A -> B) axis
        if AP_coord <= 0
            distance = norm(A.xy);
        elseif AP_coord >= AB_norm
            distance = norm(B.xy);
        else
            distance = norm(A.xy + AP_coord .* AB_unit);
        end
    end
    
    function update_info_about_D()
        if norm(D.xy) < min_geod_distances(D.id)
            min_geod_distances(D.id) = norm(D.xy);
            first_triangles(D.id, :) = first_triangle;
            first_angles(D.id) = theta;
        end
    end

    function adjacent_faces_idx = get_adjacent_faces_idx(vertex_id)
        adjacent_faces_idx = vertexAttachments(TR, vertex_id);
        adjacent_faces_idx = adjacent_faces_idx{:};
    end
end
