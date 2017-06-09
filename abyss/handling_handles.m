% Params
PARAMS = struct();

PARAMS.protocol_name = 'Waves';
PARAMS.subject_name = 'Jordy03';
PARAMS.study_name = 'jordy_spont03';
PARAMS.channels_idx = 32:99;
PARAMS.sampling_rate = 250;

PARAMS.visualize.flattening = false;

PARAMS.wave_defaults.half_width = 0.002;
PARAMS.wave_defaults.angular_half_width = pi/4;

PARAMS.max_distance = 0.02;


% load cortex
cortex = from_bst_get_surface('cortex_15002V', PARAMS);
TRI = triangulation(cortex.Faces, cortex.Vertices);

bad_edge = [775 809];

three_faces_idx = cell2mat(edgeAttachments(TRI, bad_edge));
three_faces = cortex.Faces(three_faces_idx, :);
three_faces_tri = triangulation(three_faces, cortex.Vertices);
trimesh(three_faces_tri);

some_more_faces = neighbors(TRI, three_faces_idx');

TRI_edges = edges(TRI);
edge_attachments = edgeAttachments(TRI, TRI_edges);
edge_attachment_cnts = cellfun(@length, edge_attachments);


edges_with_3 = TRI_edges(edge_attachment_cnts == 3, :);
for edge = edges_with_3'
    three_faces_idx = cell2mat(edgeAttachments(TRI, edge'));
    disp(neighbors(TRI, three_faces_idx'));
end

y = sort(edge_attachment_cnts(:));
p = find([true;diff(y)~=0;true]);
values = y(p(1:end-1));
instances = diff(p);
disp(values');
disp(instances');


cortex_smooth = from_bst_get_surface(...
    'cortex_15002V 50pct smoothed', PARAMS);
TRI_smooth = triangulation(cortex_smooth.Faces, cortex_smooth.Vertices);
vertices_of_edges_w_4_smooth = cortex_smooth.Vertices(edges_w_4_neighbours(:), :);
hold on;
trimesh(TRI_smooth);
scatter3(...
    vertices_of_edges_w_4_smooth(:, 1), ...
    vertices_of_edges_w_4_smooth(:, 2), ...
    vertices_of_edges_w_4_smooth(:, 3), 100, 'r', 'filled');
hold off;


edges_w_4_neighbours_idx = find(edge_attachment_cnts == 4);
edges_w_4_neighbours = TRI_edges(edges_w_4_neighbours_idx, :);

an_edge = edges_w_4_neighbours(1, :);
four_faces_idx = cell2mat(edgeAttachments(TRI, an_edge));
four_faces = cortex.Faces(four_faces_idx, :);
four_faces_tri = triangulation(four_faces, cortex.Vertices);
trimesh(four_faces_tri, 'FaceColor', 'r');

trimesh(four_faces, cortex.Vertices(:, 1), cortex.Vertices(:, 2), ...
    cortex.Vertices(:, 3), 'FaceColor', 'r');

four_faces_edges = vertcat(...
    four_faces(:, [1, 2]), four_faces(:, [2, 3]), four_faces(:, [3, 1]));
more_faces_idx = edgeAttachments(TRI, four_faces_edges);
more_faces_idx = setdiff(unique([more_faces_idx{:}]), four_faces_idx);
more_faces = cortex.Faces(more_faces_idx, :);
more_faces_tri = triangulation(more_faces, cortex.Vertices);
trimesh(more_faces_tri);





TRI_high = triangulation(cortex_high.Faces, cortex_high.Vertices);
TRI_high_edges = edges(TRI_high);
high_edge_attachments = edgeAttachments(TRI_high, TRI_high_edges);
high_edge_attachment_cnts = cellfun(@length, high_edge_attachments);

y = sort(high_edge_attachment_cnts(:));
p = find([true;diff(y)~=0;true]);
values = y(p(1:end-1));
instances = diff(p);
disp(values');
disp(instances');


cortex_50004V
















