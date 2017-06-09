function [faces_idx, vertices_idx] = find_artifacts(cortex)

    faces_idx = []; % faces to remove
    vertices_idx = []; % vertices to remove

    TR = triangulation(cortex.Faces, cortex.Vertices);
    TR_edges = edges(TR);
    edge_attachments = edgeAttachments(TR, TR_edges);
    edge_attachment_cnts = cellfun(@length, edge_attachments);
    
    
    % Remove dangling triangles
    edges_with_3 = TR_edges(edge_attachment_cnts == 3, :);
    for edge = edges_with_3'
        three_faces_idx = cell2mat(edgeAttachments(TR, edge'));
        adjacent_faces = neighbors(TR, three_faces_idx');
        dangling_face_id = three_faces_idx(sum(isnan(adjacent_faces), 2) == 2);
        dangling_vertex_id = setdiff(cortex.Faces(dangling_face_id, :), edge);
        
        faces_idx = [faces_idx dangling_face_id];
        vertices_idx = [vertices_idx dangling_vertex_id];
    end
    
    % Remove dangling tetrahedra
    edges_with_4 = TR_edges(edge_attachment_cnts == 4, :);
    for edge = edges_with_4'
        four_faces_idx = cell2mat(edgeAttachments(TR, edge'));
        four_faces = cortex.Faces(four_faces_idx, :);
        other_vertices = four_faces(~ismember(four_faces, edge));
        
        for i = 1:3 
        for j = (i+1):4
            if isConnected(TR, other_vertices([i, j])');
                % There might be a dangling triangle dangling from a
                % dangling tetrahedron (love this sentence). So when we
                % check that an edge connecting two out of four faces is
                % not attached to any faces except for those of the
                % potential tetrahedron we also do not mind if it is
                % attached to a dangling triangle found above.
                other_two_faces_idx = setdiff(cell2mat( ...
                    edgeAttachments(TR, other_vertices([i, j])')), faces_idx);
                if isempty(setdiff(...
                        cortex.Faces(other_two_faces_idx, :), ...
                        [edge' other_vertices([i, j])' faces_idx]))
                    tetra_faces_idx = vertexAttachments(TR, other_vertices([i, j]));
                    tetra_faces_idx = [tetra_faces_idx{:}];
                    faces_idx = [faces_idx tetra_faces_idx];
                    vertices_idx = [vertices_idx  other_vertices([i, j])'];                    
                end
            end
        end
        end
    end
    
    % Check that we found all the artifacts
    Faces = cortex.Faces;
    Faces(faces_idx, :) = [];
    TR_new = triangulation(Faces, cortex.Vertices);
    TR_new_edges = edges(TR_new);
    edge_attachments_new = edgeAttachments(TR_new, TR_new_edges);
    edge_attachment_new_cnts = cellfun(@length, edge_attachments_new);
    
    y = sort(edge_attachment_new_cnts(:));
    p = find([true;diff(y)~=0;true]);
    values = y(p(1:end-1));
    instances = diff(p);
    % disp([values instances]);
    
    assert(all(ismember(values, [1 2])), 'Could not find all the artifacts. Aborting.');
    
%     faces_idx{1} = cell2mat(edgeAttachments(TR_new, edge));
%     for i = 2:5
%         previous_layer_vertices = unique(Faces(faces_idx{i-1}, :));
%         their_adjacent_faces = vertexAttachments(TR_new, ...
%             previous_layer_vertices);
%         faces_idx{i} = setdiff([their_adjacent_faces{:}], [faces_idx{1:i-1}]);
%     end
%     
%     trimesh(Faces(faces_idx{1}, :), ...
%             cortex.Vertices(:, 1), ...
%             cortex.Vertices(:, 2), ...
%             cortex.Vertices(:, 3), ...
%             'FaceColor', 'g');
%         
%     trimesh(Faces(faces_idx{2}, :), ...
%             cortex.Vertices(:, 1), ...
%             cortex.Vertices(:, 2), ...
%             cortex.Vertices(:, 3), ...
%             'FaceColor', 'r');
%     
%     trimesh(cortex.Faces(other_two_faces_idx, :), ...
%             cortex.Vertices(:, 1), ...
%             cortex.Vertices(:, 2), ...
%             cortex.Vertices(:, 3), ...
%             'FaceColor', 'r');
end