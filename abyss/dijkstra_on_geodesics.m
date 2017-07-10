function [dist, second_on_path] = dijkstra_on_geodesics(start_vertex_id, ...
    geodesic_distances)
    % Use Dijkstra on geodesics to find shortest paths

    vertex_cnt = length(geodesic_distances);
    % dict with current shortest distance from the central vertex
    dist = inf(length(geodesic_distances), 1);
    dist(start_vertex_id) = 0;
    % list of vertices on which Dijkstra hasn't spread yet
    unvisited = 1:vertex_cnt;
    % dict containing second vertex on the shortest path
    second_on_path = -1*ones(vertex_cnt, 1);
    vertices_reachable_by_geodesics = ...
        cell2mat(geodesic_distances{start_vertex_id}.keys);
    second_on_path(vertices_reachable_by_geodesics) = ...
        vertices_reachable_by_geodesics;

    while ~isempty(unvisited)
        [~, unvisited_id] = min(dist(unvisited));
        vertex_id = unvisited(unvisited_id);
        unvisited(unvisited_id) = [];
        edges = geodesic_distances{vertex_id}.keys;
        edge_cnt = length(edges);
        for i = 1:edge_cnt
            alt_dist = ...
                dist(vertex_id) + geodesic_distances{vertex_id}(edges{i});
            if alt_dist < dist(edges{i})
                dist(edges{i}) = alt_dist;
                if vertex_id ~= start_vertex_id
                    second_on_path(edges{i}) = second_on_path(vertex_id); 
                end
            end
        end
    end                
end                