%% Params
PARAMS = struct();

PARAMS.protocol_name = 'Waves';
PARAMS.subject_name = 'Jordy03';
PARAMS.study_name = 'jordy_spont03';
PARAMS.channels_idx = 32:99;
PARAMS.sampling_rate = 250;


PARAMS.max_distance = 0.04;
PARAMS.forward.name = 'Overlapping spheres';



%% Load and procees the cortex and the measurements

% Load and flatten the cortical sheet or load if precalculated. Will take
% several hours if not precalculated.
PARAMS.cortex.low.comment = 'cortex_15002V';
PARAMS.visualize.flattening = false;
cortex = flatten_cortex_or_load(PARAMS, 'low');

% Load spike trials
PARAMS.first_peak.comment = 'spike_peak1';
meas = from_bst_get_measurements(PARAMS.first_peak.comment, PARAMS);

%% Generate waves and project them onto sensors
PARAMS.wave.half_width = 0.005; % (m)
PARAMS.wave.angular_half_width = pi/4;
PARAMS.wave.duration = 0.02; % (s)
PARAMS.wave.directions = [NaN pi/4*[0:7]];
PARAMS.wave.speeds = 0.1:0.1:1; % (mm/ms = m/s)

sensor_waves = generate_or_load_waves_on_sensors(cortex, PARAMS, 'low');

%% Evaluate how well each wave is able to describe our measurements
% Uses a lot of memory so close everything, esp. the browser


PARAMS.scores.intercept = false;
PARAMS.scores.equalize_timewise = true;
scores = evaluate_waves(sensor_waves, meas, PARAMS);
intercept = PARAMS.scores.intercept;
equalize_timewise = PARAMS.scores.equalize_timewise;
save('cortex_15002V_Jordy03_15K_4cm_scores_1', ...
    'intercept', 'equalize_timewise', 'PARAMS', 'scores', '-v7.3');

PARAMS.scores.intercept = true;
PARAMS.scores.equalize_timewise = true;
scores = evaluate_waves(sensor_waves, meas, PARAMS);
intercept = PARAMS.scores.intercept;
equalize_timewise = PARAMS.scores.equalize_timewise;
save('cortex_15002V_Jordy03_15K_4cm_scores_2', ...
    'intercept', 'equalize_timewise', 'PARAMS', 'scores', '-v7.3');

scores = load('cortex_15002V_Jordy03_15K_4cm_scores_1');
PARAMS.scores.intercept = scores.intercept;
PARAMS.scores.equalize_timewise = scores.equalize_timewise;
scores = scores.scores;

%% Visualize scores in Brainstorm
PARAMS.scores.visualize.collapse_start_time = false;
PARAMS.scores.visualize.collapse_meas = false;
move_scores_to_brainstorm(scores, PARAMS);

PARAMS.scores.visualize.collapse_start_time = true;
PARAMS.scores.visualize.collapse_meas = true;
PARAMS.scores.visualize.start_time_inds = [21:26];
move_scores_to_brainstorm(scores, PARAMS);

%% Visually find similar trials
% Use second kind of maps from above (one map per trial)
trial_ids = [1 4 8];

%% Find wave parameters common to the selected trials
% Take best (vertex, speed, direction), see which overlap
PARAMS.candidate_waves.start_time_inds = [21:26];
PARAMS.candidate_waves.percentile = 99;
[candidate_wave_inds, candidate_wave_start_times] = ...
    select_candidate_waves(scores, trial_ids, PARAMS);




%% Do the lasso shake
PARAMS.lasso.equalize_timewise = false;
PARAMS.lasso.CV_plot = false;
PARAMS.lasso.DFmax = Inf;
PARAMS.lasso.k_for_k_fold = 10;

waves_lassoed = arrayfun(@(i) ...
    lasso_the_waves(sensor_waves, candidate_wave_inds, ...
        candidate_wave_start_times, meas(trial_ids), PARAMS), ...
    [1:10], 'un', 0);

common_lassoed_waves = cellfun(@(w) ...
    find(prod(cell2mat(w) > 0, 2)), waves_lassoed, 'un', 0);
unique_wave_ids = unique(cell2mat(common_lassoed_waves')');
unique_wave_cnts = arrayfun(@(ind) ...
    nnz(cell2mat(common_lassoed_waves')' == ind), unique_wave_ids);

common_candidate_inds = unique_wave_ids(unique_wave_cnts >= 8);
final_wave_inds = candidate_wave_inds(common_candidate_inds);
final_wave_start_times = ...
    candidate_wave_start_times(common_candidate_inds, :);

final_wave_coeffs = fit_multiple_waves(sensor_waves, ...
    final_wave_inds, final_wave_start_times, meas(trial_ids), PARAMS);

% Take linear combinations, send to Brainstorm
[final_waves, final_time_points] = combine_waves(...
    final_wave_coeffs, final_wave_inds, final_wave_start_times, ...
    cortex, PARAMS);

for i = 1:length(trial_ids)
    comment = sprintf('Final solution for trial %d', trial_ids(i));
    to_bst_wave(final_waves{i}, comment, final_time_points{i}, PARAMS);
end

% Now, only directed waves
[final_waves2, final_time_points2] = combine_waves(...
    cellfun(@(x) x(d~=1), final_wave_coeffs, 'un', 0), ...
    final_wave_inds(d~=1), ...
    final_wave_start_times(d~=1, :), ...
    cortex, PARAMS);

for i = 1:length(trial_ids)
    comment = sprintf('Final solution for trial %d - directed part', trial_ids(i));
    to_bst_wave(final_waves2{i}, comment, final_time_points2{i}, PARAMS);
end

[v, s, d] = ind2sub(...
    [size(sensor_waves, 1) size(sensor_waves, 2), size(sensor_waves, 3)], ...
    final_wave_inds);

%% Playing with python
PARAMS.python.executable = 'E:\miniconda-windows\envs\matlab\python.exe';







