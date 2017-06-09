function [min_geod_distance, first_triangles, first_angles] = ...
    LOS_for_one_vertex(mesh, S_id)
% Calculates minimal-geodesic paths from the vertex with index v_id not 
% longer than max_distance.
% Also retains information about the first triangle that the path goes
% through in the form (for one destination vertex) [S_id, A_id, B_id] and
% angles between SA and the first part of the path.

    % Plot for debugging
    % DEBUG: fig = figure; hold on;

    % Array allocation
    vertex_cnt = length(mesh.Vertices);
    min_geod_distance = inf(vertex_cnt, 1);
    min_geod_distance(S_id) = 0;
    first_triangles = cell(vertex_cnt, 1);
    first_angles = NaN(vertex_cnt, 1);
    first_angles(S_id) = 0; 

    % Initialze some of the variables
    vertex_template = struct('id', NaN, 'xyz', NaN, 'xy', NaN);
    [A, B, C, D, S] = deal(vertex_template);
    S.id = S_id; S.xyz = mesh.Vertices(S_id, :); S.xy = [0, 0];
    stack = {};
    last_element_id = 0;
    TRI = triangulation(mesh.Faces, mesh.Vertices);

    % The main loop
    adjacent_triangles = mesh.Faces(any(mesh.Faces == S_id, 2), :);
    for Tau = adjacent_triangles'
        [first_triangle, theta_max] = process_first_triangle();
        theta_min = 0;
        push_to_stack(S, A, B, theta_min, theta_max);
        while ~isempty(stack)
            [C, A, B, theta_min, theta_max] = pop_from_stack();
            % DEBUG: disp(num2str([C.id A.id B.id]));
            [D, theta] = extend_over_AB();
            if isnan(D.id)
                continue;
            end
            if theta <= theta_max
                push_to_stack(A, D, B, max(theta, theta_min), theta_max);
            end
            if theta >= theta_min
                push_to_stack(B, A, D, theta_min, min(theta, theta_max));
            end
            if theta <= theta_max && theta >= theta_min
                update_info_about_D();
            end
        end
    end

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
        min_geod_distance(A.id) = norm(A.xyz - S.xyz);
        min_geod_distance(B.id) = norm(B.xyz - S.xyz);
        first_triangles{A.id} = first_triangle;
        first_triangles{B.id} = first_triangle;
        first_angles(A.id) = 0;
        first_angles(B.id) = theta_max;
    end

    function C_xy = map_triangle_to_the_plane(A, B, C)
    % Returns planar coordianates for the point C such that: 1) triangles
    % with planar and with spatial coordinates of A, B and C are congruent;
    % 2) C.xy is on the right of A.xy -> B.xy

        ABC = angle_between_3d_vectors(A.xyz - B.xyz, C.xyz - B.xyz);
        R = [cos(ABC), -sin(ABC); sin(ABC), cos(ABC)]; % Rotation matrix
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
    
    function push_to_stack(C, A, B, theta_min, theta_max)
        stack{end+1} = {C, A, B, theta_min, theta_max}; 
        % DEBUG: triplot([1 2 3], [C.xy(1), A.xy(1), B.xy(1)]', [C.xy(2), A.xy(2), B.xy(2)]');
        
    end

    function [C, A, B, theta_min, theta_max] = pop_from_stack()
        [C, A, B, theta_min, theta_max] = stack{end}{:};
        stack(end) = [];
    end

    function [D, theta] = extend_over_AB()
        ABC_and_ABD = edgeAttachments(TRI, A.id, B.id);
        if length(ABC_and_ABD{:}) == 2
            % D.id = setdiff(mesh.Faces(ABC_and_ABD{1}, :), [A.id, B.id, C.id]);
            
            vertices_id = mesh.Faces(ABC_and_ABD{1}, :);
            ABC_vert_id = [A.id, B.id, C.id];
            D.id = vertices_id(~ismembc(vertices_id(:), sort(ABC_vert_id(:))));
            
            D.xyz = mesh.Vertices(D.id, :);
            D.xy = map_triangle_to_the_plane(A, B, D);
            theta = atan2(D.xy(2), D.xy(1));
        elseif length(ABC_and_ABD{:}) == 1
            D = vertex_template;
            theta = NaN;
        elseif length(ABC_and_ABD{:}) > 2
            error('Handles in the mesh. Damn you brainstorm!');
        end
           
    end
    
    function update_info_about_D()
        if norm(D.xyz - S.xyz) < min_geod_distance(D.id)
            min_geod_distance(D.id) = norm(D.xyz - S.xyz);
            first_triangles{D.id} = first_triangle;
            first_angles(D.id) = theta;
        end
    end
end
