function angles = calculate_angles_of_paths(cortex, As, Bs, first_angles)
    
    % First, for each vertex find edges of adjacent faces that are opposite
    % the said vertex.
    TR = triangulation(cortex.Faces, cortex.Vertices);
    vertex_cnt = size(cortex.Vertices, 1);
    adj_face_ids = vertexAttachments(TR);
    adj_faces = cellfun(@(f_ids) cortex.Faces(f_ids, :)', adj_face_ids, 'un', 0);
    adj_edges = cellfun(@(f, i) reshape(f(f ~= i), 2, size(f, 2))', ...
        adj_faces, num2cell(1:vertex_cnt)', 'un', 0);
    
    % Now, sort the edges, so that [k, l] goes before [l, k]
    for v_id = 1:vertex_cnt
        e = adj_edges{v_id};
        for i = 2:size(e, 1)
            % Find the edge that contains the end of the edge i-1
            [r, c] = find(e == e(i-1, 2));
            c = c(r ~= i-1); r = r(r ~= i-1);
            % Exchange rows i and r and switch columns if needed so that
            % the end of edge i-1 is the start of the new edge i
            if c == 1, c = [1, 2]; else c = [2, 1]; end
            e([i, r], :) = e([r, i], c);
        end
        adj_edges{v_id} = e;
    end
    
    % For each vertex i and edge [k, l] find the angle (k, i, l)
    edge_angles = cell(vertex_cnt, 1);
    for v_id = 1:vertex_cnt
        B = cortex.Vertices(v_id, :);
        e = adj_edges{v_id};
        edge_cnt = size(e, 1);
        edge_angles{v_id} = zeros(edge_cnt, 1);
        for i = 1:edge_cnt
            A = cortex.Vertices(e(i, 1), :);
            C = cortex.Vertices(e(i, 2), :);
            a = A - B;
            b = C - B;
            edge_angles{v_id}(i) = 2 * atan( ...
                norm(a*norm(b) - norm(a)*b) ...
              / norm(a * norm(b) + norm(a) * b));
        end
    end
    
    % Now, cumsum the angles and add zero at the beginning for future
    % convenience. This will correspond to angles of starts of edges
    vertex_angles = cellfun(@(a) [0; cumsum(a)], edge_angles, 'un', 0);
    
    % Paths from i to j first go via the triangle (As(i,j), i, Bs(i,j)).
    % Let's call O the point where the path intersects with AB. Then angle
    % (As(i,j), i, O) is equal to first_angles(i,j). To assign an angle to
    % each path we need to find the angle we assinged to As(i,j) and to
    % establish whether A->B is a counterclockwise direction or not.
    [r, c, v_angles, v_ccw] = deal([]);
    As = As'; Bs = Bs';
    for v_id = 1:vertex_cnt
        [j, i, v_A] = find(As(:, v_id));
        [~, ~, v_B] = find(Bs(:, v_id));     
        [~, edge_n_A] = ismember(v_A, adj_edges{v_id}(:, 1));
        [~, edge_n_B] = ismember(v_B, adj_edges{v_id}(:, 1));
        r = [r; v_id*i];
        c = [c; j];
        v_angles = [v_angles; vertex_angles{v_id}(edge_n_A)];
        adj_edge_cnt = size(adj_edges{v_id}, 1);
        v_ccw = [v_ccw; ...
            (mod(edge_n_A + 1, adj_edge_cnt) ...
          == mod(edge_n_B, adj_edge_cnt)) ... % 1 if A comes before B
          - (mod(edge_n_A, adj_edge_cnt) ...
          == mod(edge_n_B + 1, adj_edge_cnt))]; % -1 if A comes after B
    end
    angles = sparse(r, c, v_angles, vertex_cnt, vertex_cnt);
    ccw = sparse(r, c, v_ccw, vertex_cnt, vertex_cnt);
    As = As'; Bs = Bs';
    
    % At this point: 
    % angles(i, j) is the angle we assigned to vertex A of the triangle
    %              (ABi)
    % ccw(i, j) is 1 if we assigned a larger angle to B than to A therefore
    %           we need to add the angle between the path and iA. ccw(i, j) 
    %           is -1 otherwise
    % first_angles(i, j) is that angle between the path and iA
    
    % Now we either add or subtract first_angles to/from angles(i, j)
    % depending on ccw(i, j) and then normalize by the sum of edge angles
    % around vertex i which can be found in vertex_angles{i}(end)
    angles = angles + ccw .* first_angles;
    normalizers = cellfun(@(a) a(end), vertex_angles);
    normalizers(normalizers == 0) = 1;
    angles = bsxfun(@rdivide, angles, normalizers) * 2*pi;
    angles = spfun(@(phi) mod(phi + pi, 2*pi) - pi, angles);
    
end