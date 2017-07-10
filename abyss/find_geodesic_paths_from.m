function LOS_output = find_geodesic_paths_from(vertex_id, bst_surf, max_distance)
    
    function D_xy = map_tri_to_the_plane(A, B, D, A_xy, B_xy)
        
        % Return x,y coordianate of the point D_xy such that:
        % 1) triABDs ABD and (A_xy, B_xy, D_xy) are congruent;
        % 2) D_xy is on the right of A_xy -> B_xy

        ABD = atan2(norm(cross(A-B, D-B)), dot(A-B, D-B));
        % cosABD = np.clip(dot(A-B, D-B)/norm(A-B)/norm(D-B), -1, 1)
        % ABD = arccos(cosABD)
        % R - rotation matrix
        R = [cos(ABD), -sin(ABD); sin(ABD), cos(ABD)];
        D_xy = B_xy' + R * (A_xy-B_xy)' ./ norm(A_xy-B_xy) .* norm(D-B);
        D_xy = D_xy';
    end

    function D_id = extend_over_edge(A_id, B_id, C_id)
        % Returns D_id != C_id such that (A_id, B_id, D_id) is a triangle
        % in faces. If no such D_id exists returns empty array

        [faces_with_A_id, ~] = find(bst_surf.Faces == A_id);
        [faces_with_B_id, ~] = find(bst_surf.Faces == B_id);
        [faces_with_C_id, ~] = find(bst_surf.Faces == C_id);
        new_face_row_id = setdiff( ...
            intersect(faces_with_A_id, faces_with_B_id), ...
            faces_with_C_id);
        D_id = setdiff(bst_surf.Faces(new_face_row_id, :), [A_id, B_id]);
    end

        disp(vertex_id);
        
        vertex = bst_surf.Vertices(vertex_id, :);
        [faces_idx_with_vertex, ~] = find(bst_surf.Faces == vertex_id);
        faces_with_vertex_cnt = length(faces_idx_with_vertex);
        
        geodesic_distances = containers.Map( ...
            'KeyType', 'int32', 'ValueType', 'double');
        first_faces = containers.Map( ...
            'KeyType', 'int32', 'ValueType', 'any');
        Thetas = containers.Map( ...
            'KeyType', 'int32', 'ValueType', 'double');        
        
        for Tau_id = 1:faces_with_vertex_cnt
            
            fprintf('Tau_id = %d\n', Tau_id);
            
            Tau_face_id = faces_idx_with_vertex(Tau_id);
            Tau = bst_surf.Faces(Tau_face_id, :);
            C_id =  vertex_id;
            A_and_B_id = Tau(Tau ~= C_id);
            A_id =  A_and_B_id(1);
            B_id =  A_and_B_id(2);
            A = bst_surf.Vertices(A_id, :);
            B = bst_surf.Vertices(B_id, :);
            Tau = [A_id B_id C_id];
            
            geodesic_distances(A_id) = norm(vertex-A);
            if first_faces.isKey(A_id)
                first_faces.remove(A_id);
                Thetas.remove(A_id);
            end
            
            A_xy = [norm(vertex-A), 0];
            B_xy = map_tri_to_the_plane(A, vertex, B, A_xy, [0, 0]);
            theta_min = 0;
            theta_max = atan2(B_xy(2), B_xy(1));
            
            chain_stack = {};
            
            while 1==1
                
                % D_id_arr = extend_over_edge(bst_surf.Faces, A_id, B_id, C_id);
                D_id_arr = extend_over_edge(A_id, B_id, C_id);                

                % For each possible extensions (yeah, there can be
                % multiple) update paths and add new edges that can still 
                % be seen to the stack
                for D_id_id = 1:length(D_id_arr)
                    D_id = D_id_arr(D_id_id);
                    D = bst_surf.Vertices(D_id, :);
                    if norm(D-vertex) >= max_distance
                        continue;
                    end
                    D_xy = map_tri_to_the_plane(A, B, D, A_xy, B_xy);
                    
                    theta = atan2(D_xy(2), D_xy(1));
                    if (theta_min < theta) && (theta < theta_max)
                        if ~geodesic_distances.isKey(D_id) || ...
                                (norm(D_xy) < geodesic_distances(D_id))
                            geodesic_distances(D_id) = norm(D_xy);
                            first_faces(D_id) = Tau;
                            Thetas(D_id) = theta;
                        end
                    end
                    if theta > theta_min
                        chain_stack{end+1} = {A_id, A_xy, D_id, D_xy, B_id, ...
                                        theta_min, min(theta, theta_max)};
                    end
                    if theta < theta_max
                        chain_stack{end+1} = {D_id, D_xy, B_id, B_xy, A_id, ...
                                        max(theta, theta_min), theta_max};
                    end                    
                end

                if ~isempty(chain_stack)
                    [A_id, A_xy, B_id, B_xy, C_id, ...
                        theta_min, theta_max] = chain_stack{end}{:};
                    A = bst_surf.Vertices(A_id, :);
                    B = bst_surf.Vertices(B_id, :);
                    chain_stack(end) = [];
                    continue;
                else
                    break;
                end
                
            end
        end
    LOS_output = {geodesic_distances, first_faces, Thetas};
end