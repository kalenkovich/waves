function [bst_surf_reduced, index_mapping] = reduce_surface(bst_surf, ...
    dipole_idx, max_distance)
    
    dipole_locations = num2cell(bst_surf.Vertices(dipole_idx, :), 2);
    
    % Distances to a closest dipole
    vectors_to_dipoles = cellfun( ...
        @(loc) bsxfun(@minus, bst_surf.Vertices, loc), dipole_locations, 'un', 0);
    distance_to_dipoles = cellfun( ...
        @(vectors) sqrt(sum(vectors.^2, 2)), vectors_to_dipoles, 'un', 0);
    min_dist_to_dipoles = min(cell2mat(distance_to_dipoles'), [], 2);

    % Remove vertices that are farther than  from all the dipoles
    too_far_vertices_idx = find(min_dist_to_dipoles >= max_distance);
    [bst_surf_reduced, index_mapping] = ...
        remove_vertices(bst_surf, too_far_vertices_idx);
    first_dipole_reduced_idx = index_mapping(dipole_idx(1));
    
    % Create the adjacency matrix and find which connected component the
    % first dipole belongs too
    faces = bst_surf_reduced.Faces;
    reduced_vertex_cnt = length(bst_surf_reduced.Vertices);
    adj_mat = zeros(reduced_vertex_cnt);
    for face_idx = 1:length(faces)
        aux = perms([1 2 3]);
        aux = aux(:, 1:2);
        face = faces(face_idx, :);
        for i = 1:6
            adj_mat(face(aux(i, 1)), face(aux(i, 2))) = 1;
        end
    end
    g = graph(adj_mat);
    c = conncomp(g);
    unconnected_vertices_idx = find(c ~= c(first_dipole_reduced_idx));
    
    % Remove vertices that are too far or that become disconnected once we
    % remove the too-far ones
    unconnected_vertices_original_idx = ...
        arrayfun(@(x) find(index_mapping == x, 1), unconnected_vertices_idx);
    too_far_or_disconnected_idx = ...
        union(too_far_vertices_idx, unconnected_vertices_original_idx);
    [bst_surf_reduced, index_mapping] = ...
        remove_vertices(bst_surf, too_far_or_disconnected_idx);
end

function [bst_surf_reduced, index_mapping] = remove_vertices(bst_surf, vertices_to_remove_idx)
    vertex_cnt = length(bst_surf.Vertices);
    vertices_to_remove_cnt = length(vertices_to_remove_idx);
    vertices_to_retain_idx = sort(setdiff(1:vertex_cnt, vertices_to_remove_idx));
    vertices_to_retain_cnt = vertex_cnt - vertices_to_remove_cnt;
    
    index_mapping = 1:vertex_cnt;
    index_mapping(vertices_to_retain_idx) = 1:vertices_to_retain_cnt;
    index_mapping(vertices_to_remove_idx) = NaN; 
    
    bst_surf_reduced = bst_surf;
    bst_surf_reduced.Vertices = bst_surf_reduced.Vertices(vertices_to_retain_idx, :);
    bst_surf_reduced.Faces = arrayfun(@(old_vertex_id) index_mapping(old_vertex_id), ...
        bst_surf_reduced.Faces);
    bst_surf_reduced.Faces(any(isnan(bst_surf_reduced.Faces), 2), :) = [];
end




    