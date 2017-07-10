% Load, change the names
cortex = load('for_evgenii3');
cortex.Vertices = cortex.vertices;
cortex.Faces = cortex.faces;
cortex = rmfield(cortex, {'vertices', 'faces'});

%% LOS distances
vertex_cnt = size(cortex.Vertices, 1);
% min_geod_distance = inf(vertex_cnt);
min_geod_distances = cell(vertex_cnt, 1);
first_triangles = cell(vertex_cnt, 1);
first_angles = cell(vertex_cnt, 1);
original_idx = cell(vertex_cnt, 1);

max_distance = 1000;

tic;
for v_id = 1:vertex_cnt
    [min_geod_distances{v_id}, first_triangles{v_id}, ...
     first_angles{v_id}] = ...
        LOS_for_one_vertex(cortex, v_id, ...
            max_distance);
    original_idx{v_id} = (1:vertex_cnt)';
end
toc;

%% Compare LOS distances

for i = 1:size(cortex.vertex_pairs, 1)
    vertex_pair = cortex.vertex_pairs(i, :);
    my_distance = min_geod_distances{vertex_pair(1)}(vertex_pair(2));
    his_distance = cortex.distances_LOS(i);
    if (abs(my_distance-his_distance) > 100*eps(my_distance)) ...
        || (isnan(his_distance) && my_distance ~= inf) ...
        || (~isnan(his_distance) && my_distance == inf)
        disp(vertex_pair);
        disp(cortex.distances_LOS(i));
        disp(min_geod_distances{vertex_pair(1)}(vertex_pair(2)));
    end
end

%% Floyd distances

% Prepare triplets (row index, columnd index, value) for sparse arrays with
% distances and angles
destination_cnts = cellfun(@length, min_geod_distances);
i = cell2mat(arrayfun(@(v_id, dest_cnt) ones(dest_cnt, 1) * v_id, ...
    (1:vertex_cnt)', destination_cnts, 'un', 0));
j = cell2mat(original_idx);
v_distances = cell2mat(min_geod_distances);
v_angles = cell2mat(first_angles);


inf_ids = v_distances == inf;
i(inf_ids) = [];
j(inf_ids) = [];
v_distances(inf_ids) = [];
v_angles(inf_ids) = [];

min_geod_distances = sparse(i, j, v_distances, vertex_cnt, vertex_cnt);
first_angles = sparse(i, j, v_angles, vertex_cnt, vertex_cnt);

tic;
[distances, first_stops] = apply_floyd(min_geod_distances, max_distance);
toc;

%% Compare Floyd distances
for i = 1:size(cortex.vertex_pairs, 1)
    vertex_pair = cortex.vertex_pairs(i, :);
    my_distance = distances(vertex_pair(1), vertex_pair(2));
    his_distance = cortex.distances_LOS_FLOYD(i);
    if (abs(my_distance-his_distance) > 100*eps(my_distance)) ...
        || (isnan(his_distance) && my_distance ~= inf) ...
        || (~isnan(his_distance) && my_distance == inf)
        disp(vertex_pair);
        disp(his_distance);
        disp(my_distance);
    end
end
