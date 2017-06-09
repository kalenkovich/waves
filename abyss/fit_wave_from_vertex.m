function [wave_fitted, fit_info] = fit_wave_from_vertex(cortex, vertex_id, ...
    meas, gain_matrix, PARAMS)

    % Cut the cortex around the vertex corresponding to the second peak. Reduce
    % the gain matrix accordingly.
    [planar_mesh, G, index_mapping] = cut_and_flatten_lscm(cortex, vertex_id, ...
        gain_matrix, PARAMS);

    % Find an optimal speed of propagation
    speed_fitted = fit_speed(planar_mesh, meas, G, PARAMS);
    speed_fitted = 0.5;

    % Find the right direction
    wave_fitted = zeros(length(cortex.Vertices), size(meas, 2));
    [wave_fitted(~isnan(index_mapping), :), fit_info] = ...
        fit_direction(planar_mesh, meas, G, speed_fitted, PARAMS);

    fit_info.reached_vertex_id = ...
        find(index_mapping == fit_info.reached_vertex_id);
end