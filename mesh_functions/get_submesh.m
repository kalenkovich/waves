function submesh = get_submesh(mesh, v_id, radius)
% Returns a submesh of mesh such that all the vertices in the submesh are
% not further than radius from the vertex with index v_id.
% The result may contain components unconnected to the v_id vertex. The
% component that does contain v_id may contain multiple boundaries (holes).

    % Find distances
    vectors_from_v_id = ...
        bsxfun(@minus, mesh.Vertices, mesh.Vertices(v_id, :));
    distances_from_v_id = sqrt(sum(vectors_from_v_id.^2, 2));
    
    % Filter vertices based on distances
    submesh.original_idx = find(distances_from_v_id <= radius);
    submesh.Vertices = mesh.Vertices(submesh.original_idx, :);
    
    % Create a mapping from id in mesh to id in submesh
    original_vertex_cnt = size(mesh.Vertices, 1);
    submesh_vertex_cnt = size(submesh.Vertices, 1);
    index_mapping = NaN(original_vertex_cnt, 1);
    index_mapping(submesh.original_idx) = 1:submesh_vertex_cnt;
    
    % Filter faces and use the mapping above so that faces now refer to
    % indices of vertices in submesh and not the original mesh
    submesh.Faces = index_mapping(mesh.Faces);
    submesh.Faces(any(isnan(submesh.Faces), 2), :) = [];
    submesh.index_mapping = index_mapping;
    
    % Save the index of v_id in submesh
    submesh.v_id_in_submesh = index_mapping(v_id);
end