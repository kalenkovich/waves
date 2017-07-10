function move_scores_to_brainstorm(scores, PARAMS)

    speed_dim = 2;
    direction_dim = 3;
    start_time_inds = PARAMS.scores.visualize.start_time_inds;
    maps = cellfun(@(s) ...
        squeeze(max(max(s.^2, [], speed_dim), [], direction_dim)), ...
        scores, 'un', 0);
    
    if PARAMS.scores.visualize.collapse_start_time
        start_time_dim = 2;
        maps = cellfun(@(m) max(m(:, start_time_inds), [], start_time_dim), ...
            maps, 'un', 0);
    end
    
    if PARAMS.scores.visualize.collapse_meas
        maps = {cell2mat(maps')};
    end
    
    trial_cnt = length(maps);
    for t_id = 1:trial_cnt
        comment = sprintf('%d-mm-wide wave scores, intercept = %d', ...
        PARAMS.wave.half_width*1000, PARAMS.scores.intercept);
        if ~PARAMS.scores.visualize.collapse_meas
            comment = [comment sptintf(', trial %d', t_id)];
        end
        to_bst_map(maps{t_id}, comment, PARAMS);
    end
end

function to_bst_map(map, comment, PARAMS)

    sources_struct = load('sources_stub.mat');
    sources_struct.ImageGridAmp = map;
    time_point_cnt = size(map, 2);
    if ~PARAMS.scores.visualize.collapse_start_time
        sources_struct.Time = [-time_point_cnt+1:0] / PARAMS.sampling_rate;
        if PARAMS.scores.visualize.collapse_meas
            sources_struct.Time = repmat(sources_struct.Time, ...
                1, size(map, 2)/time_point_cnt);
        end
    else
        sources_struct.Time = 1:size(map, 2);
    end
    sources_struct.Comment = comment;
    
    [~, iStudy] = bst_get('Study');
    db_add(iStudy, sources_struct);
end