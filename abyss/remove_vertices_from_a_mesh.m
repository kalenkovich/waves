function mesh = remove_vertices_from_a_mesh(mesh, vertex_ids)
% Remove vertices with indices vertex_ids from the mesh and removes any
% faces that contained them

    % Mapping from current indices to the new ones
    vertex_cnt = size(mesh.Vertices, 1);
    to_keep_cnt = vertex_cnt - length(vertex_ids);
    old_to_new_mapping = 1:vertex_cnt;
    old_to_new_mapping(vertex_ids) = NaN;
    old_to_new_mapping(~isnan(old_to_new_mapping)) = 1:to_keep_cnt;
    
    mesh.Vertices(vertex_ids, :) = [];
    mesh.Faces = old_to_new_mapping(mesh.Faces);
    mesh.Faces(any(isnan(mesh.Faces), 2), :) = [];
    
end