function [candidate_wave_inds, candidate_wave_start_times] = ...
    select_candidate_waves(scores, trial_ids, PARAMS)

    
    start_time_inds = PARAMS.candidate_waves.start_time_inds;
    start_time_dim = 4;
    trial_cnt = length(trial_ids);
    maps = cell(trial_cnt, 1); % for each (vertex, speed, direction) best score 
                               % disregarding the starting time
    maps_start_times = cell(trial_cnt, 1); % retain the best starting time
    for t_id_id = 1:trial_cnt
        t_id = trial_ids(t_id_id);
        [maps{t_id_id}, maps_start_times{t_id_id}] = ...
            max(scores{t_id}(:, :, :, start_time_inds).^2, ...
            [], start_time_dim);
    end

    percentile = PARAMS.candidate_waves.percentile;
    best = cellfun(@(m) find(m>=prctile(m(:), percentile)), maps, 'un', 0);


    common_best_inds = best{1};
    for t_id = 1:trial_cnt
        common_best_inds = intersect(common_best_inds, best{t_id});
    end

    % Add a non-directed wave with the same starting time to each wave
    scores_sz = size(scores{1}); % vertex_cnt x speed_cnt x direction_cnt x []
    % Best time is individual for each trial
    st = cellfun(@(sts) sts(common_best_inds), maps_start_times, 'un', 0);
    st = cell2mat(st');
    [v, s, d] = ind2sub(scores_sz(1:3), common_best_inds);
    v_s_d_st = unique([...
        v s d             st; ...
        v s ones(size(d)) st], 'rows');
    v_s_d = num2cell(v_s_d_st(:, 1:3), 1);
    candidate_wave_inds = sub2ind(scores_sz(1:3), v_s_d{:});
    candidate_wave_start_times = v_s_d_st(:, 4:end);

end