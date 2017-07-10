function [neighbours_idx, flat_coordinates] = ...
                                flatten_around_a_vertex(surf, vertex_id)
    tri = triangulation(surf.Faces, surf.Vertices);
    
    % Find all the points connected to the point with vertex_id and list
    % them in consecutive order
    neighbour_faces_idx = cell2mat(vertexAttachments(tri, vertex_id));
    neighbour_cnt = length(neighbour_faces_idx);
    neighbours_idx = points_from_faces_sorted(tri, ...
        neighbour_faces_idx, vertex_id);
    
    % Find angles between edges connecting the point with vertex_id to
    % pairs of its neigbours. Also, find norms of those edges.
    % angles(i) - angle between the edge to neighbour i and the next one
    % norms(i) - norm of the edge to neighbour i
    [angles, norms] = calculate_angles(tri, vertex_id, neighbours_idx);
    
    % Standardize angles so that they add up to 2*pi. Then do cumsum so
    % that angles become angles between future flattened points and the
    % x axis. We will also have to shift angles by 1 so that the first
    % vertex lies on the x axis.
    angles = angles / sum(angles) * (2*pi);
    angles = cumsum(angles);
    angles = circshift(angles, 1);
    
    % Find flattened points so that the vertex with vertex_id becomes (0,
    % 0), the norm of the edges stays the same and the angles are those
    % that we have just calculated
    flat_coordinates = bsxfun(@times, norms, [cos(angles), sin(angles)]);
    neighbours_idx(end+1) = vertex_id;
    flat_coordinates(end+1, :) = [0, 0];
end

function neighbour_points_sorted_idx = points_from_faces_sorted(tri, ...
                                        neighbour_faces_idx, vertex_id)
    neighbour_faces_mat =  tri.ConnectivityList(neighbour_faces_idx, :);
    neighbour_cnt = length(neighbour_faces_idx);
    neighbour_points_sorted_idx = zeros(neighbour_cnt, 1);
    neighbour_points_sorted_idx(1:2) = ...
        setdiff(neighbour_faces_mat(1, :), vertex_id);
    for i = 3:neighbour_cnt
        previous_point_id = neighbour_points_sorted_idx(i-1);
        faces_with_previous_mask = any(ismember(neighbour_faces_mat, ...
            previous_point_id), 2);
        preprevious_point_id = neighbour_points_sorted_idx(i-2);
        faces_with_preprevious_mask = any(ismember(neighbour_faces_mat, ...
            preprevious_point_id), 2);
        next_face = neighbour_faces_mat( ...
            faces_with_previous_mask & ~faces_with_preprevious_mask, :);
        neighbour_points_sorted_idx(i) = ...
            setdiff(next_face, [vertex_id, previous_point_id]);
    end
end

function [angles, norms] = calculate_angles(tri, vertex_id, neighbours_idx)
    neighbour_cnt = length(neighbours_idx);
    angles = zeros(neighbour_cnt, 1);
    norms = zeros(neighbour_cnt, 1);
    for i = 1:neighbour_cnt
        C = tri.Points(vertex_id, :);
        A = tri.Points(neighbours_idx(i), :);
        B = tri.Points(neighbours_idx(mod(i, neighbour_cnt) + 1), :);
        a = A - C;
        b = B - C;
        angles(i) = angle_between_vectors(a, b);
        norms(i) = norm(a);
    end
    
end

function angle = angle_between_vectors(a, b)
    angle = atan2(norm(cross(a, b)), dot(a, b));
end

% Testing on an already flat mesh
% tmp_tri = [ones(1, 6); 2:7; [3:7 2]]';
% tmp_points = [[0 cos(2*pi/6.*[0:5])]', ...
%               [0 sin(2*pi/6.*[0:5])]', ...
%               zeros(7, 1)];
% scatter(tmp_points(:,1), tmp_points(:,2), 10, 1:7);
% tmp_cortex = struct('Faces', tmp_tri, 'Vertices', tmp_points);
% [tmp_neighbours_idx, tmp_flat_coordinates] = ...
%                 flatten_around_a_vertex(tmp_cortex, 1);