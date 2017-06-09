function [planar_mesh, G, index_mapping] = cut_and_flatten(cortex, ...
    vertex_id, gain_matrix, PARAMS)
    
    max_distance = PARAMS.max_distance;
    cortex_tri = triangulation(cortex.Faces, cortex.Vertices);

    % Get the immediate neigbours and flatten the corresponding mesh by a)
    % keeping lengths of edges, b) normalizing the central vertex angles so
    % that they add up to 2pi
    [neighbours_idx, flat_coordinates] = flatten_around_a_vertex(cortex, ...
        vertex_id);
    
    % Initialize planar_mesh and index_mapping with the central vertex and
    % its immediate neigbours. We'll do the faces at the end.    
    planar_mesh = struct('Vertices', flat_coordinates, ...
                         'Faces', []);    
    original_vertex_cnt = size(cortex.Vertices, 1);
    index_mapping = nan(original_vertex_cnt, 1);
    index_mapping(neighbours_idx) = 1:length(neighbours_idx);
    
    % While at least one vertex on the boundary is within max_distance from
    % the central vertex (euclidian distance in the original cortex) add
    % another layer of vertices
    boundary_vertices_idx = setdiff(neighbours_idx, vertex_id);
    min_boundary_distance = get_min_boundary_distance();
    while min_boundary_distance < max_distance
        
        % New layer consists of neighbours of the vertices that have
        % already been included
        old_boundary_vertices_idx = boundary_vertices_idx;
        boundary_vertices_idx = get_closest_neighbours();
        min_boundary_distance = get_min_boundary_distance();
        
        % Add vertices connected to two old ones
        vertices_with_two_old = get_vertices_with_two_old();
        for v_id = 1:length(vertices_with_two_old)
            vertex = vertices_with_two_old(v_id);
            add_vertex_with_two_old(vertex);
        end
        draw_flat();
        
        % Add vertices connected to one old one. First, find those vertices
        % and group connected ones. Then, add connected groups ony by one.
        groups_with_one_old = get_groups_with_one_old();
        for g_id = 1:length(groups_with_one_old)
            
            % Each group is of format [old_vertex,
            % new_vertex_connect_with_two_old_ones_1,
            % some_vertices_connected_to_one_old_one,
            % new_vertex_connect_with_two_old_ones_2]
            group = groups_with_one_old{g_id};
            
            % Groups of length 3 do not contain any new vertices
            if length(group) > 3
                add_group_with_one_old(group);
                draw_flat();
            end
        end
        
    end
    
    
    
    
    
    
    % ------------- Helper functions ---------------------
    function min_boundary_distance = get_min_boundary_distance()
        vectors_from_central_vertex = bsxfun(@minus, ...
            cortex.Vertices(boundary_vertices_idx, :), ...
            cortex.Vertices(vertex_id, :));
        min_boundary_distance = ...f
            min(norms_of_rows(vectors_from_central_vertex));        
    end

    function boundary_vertices_idx = get_closest_neighbours()
        already_included_vertices = find(~isnan(index_mapping));
        neighbour_faces_idx = vertexAttachments(cortex_tri, ...
            already_included_vertices);
        neighbour_faces_idx = unique([neighbour_faces_idx{:}]);
        neighbour_vertex_idx = unique(cortex.Faces(neighbour_faces_idx, :));
        boundary_vertices_idx = setdiff(neighbour_vertex_idx, ...
            already_included_vertices);
    end
    
    function vertices_with_two_old = get_vertices_with_two_old()
        faces_with_two_old = find(sum(...
            ismember(cortex_tri.ConnectivityList, ...
                     old_boundary_vertices_idx), 2) == 2);
        faces_with_one_new = find(sum(...
            ismember(cortex_tri.ConnectivityList, ...
                     boundary_vertices_idx), 2) == 1);
        faces_with_two_old_one_new = intersect(faces_with_two_old, ...
            faces_with_one_new);
        vertices_with_two_old_one_new = unique(...
            cortex_tri.ConnectivityList(faces_with_two_old_one_new, :));
        vertices_with_two_old = setdiff(vertices_with_two_old_one_new, ...
            old_boundary_vertices_idx);
    end

    function add_vertex_with_two_old(vertex)
        faces_with_two_old = find(sum(...
            ismember(cortex_tri.ConnectivityList, ...
                     old_boundary_vertices_idx), 2) == 2);
        faces_with_this_vertex = find(sum(...
            ismember(cortex_tri.ConnectivityList, ...
                     vertex), 2) == 1);
        face_to_add = intersect(faces_with_two_old, faces_with_this_vertex);
        old_vertices = setdiff(cortex_tri.ConnectivityList(face_to_add, :), ...
            vertex);
        old_vertices_xy = planar_mesh.Vertices(index_mapping(old_vertices), :);
        sides = norms_of_rows(bsxfun(@minus, ...
            cortex_tri.Points(old_vertices, :), ...
            cortex_tri.Points(vertex, :)));        
        vertex_xy = finish_a_triangle(old_vertices_xy, sides);
        add_already_flattened_vertex(vertex, vertex_xy);
    end

    function groups_with_one_old = get_groups_with_one_old()
        % Each cell of the output looks like this: [common_old_neighbour,
        % vertex_with_two_old1, several_vertices_with_one_old,
        % vertex_with_two_old2]
        vert_with_two_cnt = length(vertices_with_two_old);
        groups_with_one_old = cell(vert_with_two_cnt, 1);
        for v_id = 1:vert_with_two_cnt
            old_vertex = old_boundary_vertices_idx(v_id);
            groups_with_one_old{v_id} = old_vertex;
            connected_idx = boundary_vertices_idx(arrayfun(@(v) ...
                isConnected(cortex_tri, old_vertex, v), ...
                boundary_vertices_idx));
            with_two_old_idx = intersect(connected_idx, ...
                vertices_with_two_old);
            groups_with_one_old{v_id}(end+1) = with_two_old_idx(1);
            while length(groups_with_one_old{v_id}) < ...
                    length(connected_idx) + 1
                faces_idx = cell2mat(edgeAttachments(cortex_tri, ...
                    old_vertex, groups_with_one_old{v_id}(end)));
                next_one_id = setdiff(...
                    intersect(...
                        unique(cortex_tri.ConnectivityList(faces_idx, :)), ...
                        boundary_vertices_idx), ...
                    groups_with_one_old{v_id});
                groups_with_one_old{v_id}(end+1) = next_one_id;
            end
        end
    end

    function add_group_with_one_old(group)
        % Calculate and normalize angles, calculate distances
        vectors = num2cell(bsxfun(@minus, ...
            cortex_tri.Points(group(2:end), :), ...
            cortex_tri.Points(group(1), :)), 2);
        angles = cellfun(@(v1, v2) angle_between_vectors(v1, v2), ...
            vectors(1:end-1, :), vectors(2:end, :));
        end_vertex_xy = planar_mesh.Vertices(...
            index_mapping([group(2) group(end)]), :);
        flat_angle = angle_between_vectors(...
            end_vertex_xy(1, :), end_vertex_xy(2, :));
        angles = cumsum(angles * flat_angle / sum(angles));
        distances = norms_of_rows(cell2mat(vectors(2:end)));
        
        % Calculate flattened coordinates
        old_vertex_xy = planar_mesh.Vertices(index_mapping(group(1)), :);
        for v_id = 1:(length(angles) - 1)
            vertex_xy = insert_into_an_angle(old_vertex_xy, end_vertex_xy, ...
                angles(v_id), distances(v_id));
            % v_id starts with one, first new vertex in group is the third
            % one so we have to do v_id+2
            add_already_flattened_vertex(group(v_id+2), vertex_xy); 
        end
    end

    function D_xy = insert_into_an_angle(C_xy, AB_xy, ...
                theta, r)
    % Find coordinates of such a point D that angle ACD is alpha, D is
    % inside the angle ACB, CD = r
        % https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula
        vectors = bsxfun(@minus, AB_xy, C_xy);
        a = [vectors(1, :) 0]; b = [vectors(2, :) 0];
        k = cross(a, b)/norm(cross(a,b));
        v = a/norm(a) * r;
        v_rot = v * cos(theta) + cross(k, v) * sin(theta) ...
            + k * dot(k, v) * (1 - cos(theta));
        D_xy = C_xy + v_rot(1:2);
        
    end

    function angle = angle_between_vectors(a, b)
        if length(a) == 2
            a(end+1) = 0; b(end+1) = 0;
            angle = angle_between_vectors(a, b);
        else
            angle = atan2(norm(cross(a, b)), dot(a, b));
        end
    end

    function add_already_flattened_vertex(vertex, vertex_xy)
        planar_mesh.Vertices(end+1, :) = vertex_xy;
        index_mapping(vertex) = length(planar_mesh.Vertices);
    end

    function draw_flat()
        figure;
        fill_faces();
        triplot(triangulation(planar_mesh.Faces, planar_mesh.Vertices));
        waitforbuttonpress;
    end

    function fill_faces()
        planar_mesh.Faces = arrayfun(...
            @(old_vertex_id) index_mapping(old_vertex_id), ...
            cortex_tri.ConnectivityList);
        planar_mesh.Faces(any(isnan(planar_mesh.Faces), 2), :) = [];
    end

    % ----------------- Visualizations -------------------
    if isfield(PARAMS, 'visualize') ...
            && isfield(PARAMS.visualize, 'flattening') ...
            && PARAMS.visualize.flattening == true
        
        % cortex and reduced
        figure; hold on;
        colors = ones(1, size(cortex.Vertices, 1));
        cortex_patch = trimesh(cortex.Faces, cortex.Vertices(:, 1), ... 
            cortex.Vertices(:, 2), cortex.Vertices(:, 3), colors);

        colors_reduced = 5*ones(1, sum(~isnan(index_mapping)));
        reduced_patch = trimesh(cortex_reduced.Faces, ... 
            cortex.Vertices(~isnan(index_mapping), 1), ...
            cortex.Vertices(~isnan(index_mapping), 2), ...
            cortex.Vertices(~isnan(index_mapping), 3), ...
            colors_reduced);
        
        hold off;
        
        % reduced with central vertex and its neigbours
        figure; hold on;
        
        colors_reduced = 2*ones(1, sum(~isnan(index_mapping)));
        reduced_patch = trimesh(cortex_reduced.Faces, ... 
            cortex.Vertices(~isnan(index_mapping), 1), ...
            cortex.Vertices(~isnan(index_mapping), 2), ...
            cortex.Vertices(~isnan(index_mapping), 3), ...
            colors_reduced);
        
        scatter3(...
            cortex_reduced.Vertices(neighbours_idx, 1), ...
            cortex_reduced.Vertices(neighbours_idx, 2), ...
            cortex_reduced.Vertices(neighbours_idx, 3), ...
            'filled');
        
        hold off;
        
        % flattened with central vertex and its neigbours
        figure; hold on;
        
        triplot(planar_mesh.Faces, ...
            planar_mesh.Vertices(:, 1), ...
            planar_mesh.Vertices(:, 2), 'm');
        scatter(flat_coordinates(:,1), flat_coordinates(:,2), 'filled');
        
        hold off;
   end
    
end

function vertex_xy = finish_a_triangle(vertices_xy, sides)
% Given coordinates of A, B and sides a, b return coordinates of C.
% There are two possible answers. Returns the one furthe from zero.
    % http://math.stackexchange.com/a/1367732
    r1 = sides(1); r2 = sides(2);
    O = mean(vertices_xy); % middle of AB
    AB = vertices_xy(2, :) - vertices_xy(1, :);
    R = norm(AB); % length of c
    av = AB./R; % unit vector along AB
    bv = [AB(2), -AB(1)]./R; % its complement
    a = (r1^2 - r2^2)/(2*R);
    b = sqrt(...
        (r1^2 + r2^2)/2 - a^2 - R^2/4);
    xy1 = O + a.*av + b.*bv;
    xy2 = O + a.*av - b.*bv;

    % Pick one that is further from zero
    if norm(xy1) > norm(xy2)
        vertex_xy = xy1;
    else
        vertex_xy = xy2;
    end
end