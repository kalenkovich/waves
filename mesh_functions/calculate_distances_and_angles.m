function cortex = calculate_distances_and_angles(cortex, PARAMS)
% For each vertex calculates lengths of the shortest paths and "angles" of
% these paths. Angles are calculated in the following way: each path first
% goes through an adjacent triangle or an edge. We assign one of the edges
% to have angle 0. Then we choose a direction and calculate an angle for a
% path as sum of the angles between the path and the 0-edge and then (if
% applicable) the angle inside the triangle that the path starts with. Such
% angles won't add up to 2pi so we then normalize them.
%
% First, for each vertex calculates the minimal-geodesic distances using
% the LOS algorithm
% (https://www.cis.upenn.edu/~cis610/geodesics-PAMI09-Schwartz.pdf)
%
% Then, using the Floyd algorithm on the results we get the shortes paths
% (also see the article cited above)
%
% While we were doing all of that we also saved enough information to
% recover the actual paths. For each vertex we can now find first legs of
% paths from it and the angle they make inside the first triangle on their
% way. Combining this information we can calculate angles of path going
% from any one vertex and normalize them so that they add up to 2pi.

    max_distance = PARAMS.max_distance;
    
    % We need to remove certain faces so that the valence of all edges is 2
    [artifacts.faces_idx, artifacts.vertices_idx] = find_artifacts(cortex);
    Faces = cortex.Faces;
    cortex.Faces(artifacts.faces_idx, :) = [];
    warning('off', 'MATLAB:triangulation:PtsNotInTriWarnId');

    %% Calculate minimal-geodesic distances
    
    % Preallocation
    vertex_cnt = size(cortex.Vertices, 1);
    min_geod_distances = cell(vertex_cnt, 1);
    first_triangles = cell(vertex_cnt, 1);
    first_angles = cell(vertex_cnt, 1);
    original_idx = cell(vertex_cnt, 1);
    
    tic;
    for v_id = 1:vertex_cnt
        submesh = get_submesh(cortex, v_id, max_distance);
        [min_geod_distances{v_id}, first_triangles{v_id}, ...
         first_angles{v_id}] = ...
            LOS_for_one_vertex(submesh, submesh.v_id_in_submesh, ...
                max_distance);
        original_idx{v_id} = submesh.original_idx;
        
        % Status update
        if mod(v_id, 100) == 0 || v_id == 1 || v_id == vertex_cnt
            fprintf('LOS: %d iterations out of %d. ', ...
                v_id, vertex_cnt);
            toc;
        end        
    end
    
    %% Clean up    
    
    % Prepare triplets for sparse arrays with distances and angles
    destination_cnts = cellfun(@length, min_geod_distances);
    i = cell2mat(arrayfun(@(v_id, dest_cnt) ones(dest_cnt, 1) * v_id, ...
        (1:vertex_cnt)', destination_cnts, 'un', 0));
    j = cell2mat(original_idx);
    v_distances = cell2mat(min_geod_distances);
    v_angles = cell2mat(first_angles);
    
    % Remove elements corresponding to unreachable pairs
    inf_ids = v_distances == inf;
    i(inf_ids) = [];
    j(inf_ids) = [];
    v_distances(inf_ids) = [];
    v_angles(inf_ids) = [];
    
    % first_triangles currently contains triplets of vertex indices in the
    % submeshes. Let's convert this to the original indices.
    non_zero_rows = @(mat) mat(any(mat, 2), :);
    first_tri_original = cellfun(@(orig, first) orig(non_zero_rows(first)), ...
        original_idx, first_triangles, 'un', 0);
    first_tri_original = cell2mat(first_tri_original);
    
    % i and j are the same, first vertex is always the same as i. What we
    % need are the second (A) and the third vertex (B)
    v_A = first_tri_original(:, 2);
    v_B = first_tri_original(:, 3);
    
    
    min_geod_distances = sparse(i, j, v_distances, vertex_cnt, vertex_cnt);
    first_angles = sparse(i, j, v_angles, vertex_cnt, vertex_cnt);
    As = sparse(i, j, v_A, vertex_cnt, vertex_cnt);
    Bs = sparse(i, j, v_B, vertex_cnt, vertex_cnt);
    
    clearvars -except min_geod_distances first_angles As Bs ...
        cortex PARAMS max_distance vertex_cnt Faces;
    
    %% Assign angles to geodesic paths
    fprintf('Assigning angles to all geodesic paths found by LOS\n');
    angles_of_geodesic_paths = calculate_angles_of_paths(cortex, ...
        As, Bs, first_angles);    
    
    %% Calculate shortest-path distances
    
    % Apply Floyd-Warshall algorithm on the graph of shortest geodesic
    % distances. Save the second end of the first leg of the shortest paths
    % so that it is possible to reconstruct the paths themselves.
    [distances, first_stops] = apply_floyd(min_geod_distances, max_distance);
    
    % It should not matter much but any shortest path with length >
    % max_distance is most likely wrong. There should not actually be any.
    distances(distances > max_distance) = 0;

    %% Assign angles to shortest paths
    [i, j, v] = find(first_stops);
    angle_ids = sub2ind(size(first_stops), i, v);
    v_angles = angles_of_geodesic_paths(angle_ids);
    angles = sparse(i, j, v_angles, vertex_cnt, vertex_cnt);
    
    %% Add results to cortex
    cortex.distances = distances;
    cortex.angles = angles;
    
    cortex.geodesic.distances = min_geod_distances;
    cortex.geodesic.first_angles = first_angles;
    cortex.geodesic.As = As;
    cortex.geodesic.Bs = Bs;
    cortex.Faces = Faces;
    cortex.PARAMS = PARAMS;
    
end
