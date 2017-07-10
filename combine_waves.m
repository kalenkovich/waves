function [combined_waves, time_points] = ...
    combine_waves(coeffs, inds, start_times, cortex, PARAMS)
    
    combined_waves = cell(size(coeffs));
    time_points = cell(size(coeffs));
    time_point_cnt = get_time_point_cnt(PARAMS);
    sampling_rate = PARAMS.sampling_rate;
    
    vertex_cnt = size(cortex.Vertices, 1);
    speed_cnt = length(PARAMS.wave.speeds);
    direction_cnt = length(PARAMS.wave.directions);
    [v, s, d] = ind2sub([vertex_cnt speed_cnt direction_cnt], inds);
    cortex_waves = arrayfun(@(v_id, s_id, d_id) ...
        generate_wave(v_id, s_id, d_id, cortex, PARAMS), v, s, d, 'un', 0);
    
    trial_cnt = length(coeffs);
    for trial_id = 1:trial_cnt
        st = start_times(:, trial_id);
        time_points{trial_id} = [min(st)-time_point_cnt : max(st)-1] ...
            / sampling_rate;
        
        qs = coeffs{trial_id};
        waves = cellfun(@(cortex_wave) ...
            padarray(cortex_wave, [0 max(st)-min(st)], 0, 'post'), ...
            cortex_waves, 'un', 0);
        waves = arrayfun(@(wave, st_id) ...
            circshift(wave{1}, st_id - min(st), 2), ...
            waves, st, 'un', 0);
        waves = arrayfun(@(wave, q) wave{1}*q, waves, qs, 'un', 0);
        combined_waves{trial_id} = sum(cat(3, waves{:}), 3);
    end
end