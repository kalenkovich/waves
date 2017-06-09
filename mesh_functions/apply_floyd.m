function [D, first_stops] = sparse_Floyd_Warshall(costs, max_distance)
% Calculates distances using Floyd-Warshall algorithm. The input costs is a
% sparse array of distances between connected vertices. first_stops(i, j)
% is a vertex that is the first_stop on the shortest path between i and j.
% It is assumed that finite costs are not higher than max_distance and that
% we are interested only in paths that are not longer than max_distance.

% Note: in sparse distant array zeros on the diagonals are actual zeros and
% zeros elsewhere stand for infinity.
    
    D = costs;
    vertex_cnt = size(D, 1);
    [r, c] = find(D);
    first_stops = sparse(r, c, c, vertex_cnt, vertex_cnt);

    tic;
    for k = 1:vertex_cnt
        % Construct D_ikj - matrix of lengths of paths i->k->j with zeros
        % on diagonal as we certainly don't need i->k->i paths.
        reachable_ids = find(D(:, k));
        reachable_dists = nonzeros(D(:, k));
        reachable_cnt = length(reachable_ids);
        i = repmat(reachable_ids, reachable_cnt, 1);
        j = repmat(reachable_ids', reachable_cnt, 1);
        j = j(:);
        v_ik = repmat(reachable_dists, reachable_cnt, 1);
        v_kj = repmat(reachable_dists', reachable_cnt, 1);
        v_kj = v_kj(:);
        D_ikj = sparse(i, j, v_ik + v_kj, vertex_cnt, vertex_cnt);
        D_ikj = D_ikj - diag(diag(D_ikj));
        
        % Update D where it makes sense to go through k. This happens in
        % two cases:
        % 1) i->k->j is shorter than finite i->j
        % 2) i->k->j is less than max_distance and i->k is infinite
        to_update_1 = spones(D_ikj) .* (D_ikj < D);
        [r, c, v] = find(D_ikj);
        r(v > max_distance) = [];
        c(v > max_distance) = [];
        v(v > max_distance) = [];
        to_update_2_1 = sparse(r, c, v, vertex_cnt, vertex_cnt);
        to_update_2_2 = D_ikj - D_ikj.*spones(D);
        to_update_2 = to_update_2_1.*to_update_2_2;
        to_update = spones(to_update_1 + to_update_2);
        D = D - D.*to_update + D_ikj.*to_update;
        first_stops = first_stops - first_stops.*to_update ...
            + to_update*diag(first_stops(:, k));
        
        if mod(k, 100) == 0 || k == 1 || k == vertex_cnt
            fprintf('Floyd-Warshall: %d iterations out of %d. ', ...
                k, vertex_cnt);
            toc;
        end
    end
    
    
end