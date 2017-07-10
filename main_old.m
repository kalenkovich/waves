%% Params
PARAMS = struct();
PARAMS.subject_name = 'Jordy03';
PARAMS.channels_idx = 32:99;
PARAMS.sampling_rate = 250;

%% House-keeping stuff with surfaces, indexes, dipole vertexes, etc.

% @var dipoles - dipole from bst
% @var bst_surf - cortex from bst used for the head model
clearvars -except dipoles bst_surf;

dipole_cnt = length(dipoles.Dipole);
dipole_idx = zeros(dipole_cnt, 1); % indices of dipoles in the whole cortex
for i = 1:dipole_cnt
    [~, dipole_idx(i)] = ...
        ismember(dipoles.Dipole(i).Loc', bst_surf.Vertices, 'rows');
end


% Switch bst_surf for a smoothed version of the cortex!!!
% @var bst_surf - now, smoothed version of the cortex. Need it because it
% does not contain handles

clearvars -except dipoles dipole_idx bst_surf dipole_cnt;

max_distance = 0.02;

% Remove vertices that are farther than max_distance from all the dipoles
% and are not connected to the first dipole. Retain mapping between
% vertices of the original and the reduced surface. 

% @var bst_surf_reduced - same format as the bst_surf, but fields Vertices
% and Faces are change to inculed only the vertices that are not too far
% @var index_mapping - an array to retain mapping between the whole cortex
% surface and its reduced version. index_mapping[n] = k if vertex number n
% in the original surface is the vertex number k in the reduced surface.
% NaN if the vertex is not in the reduced surface

[bst_surf_reduced, index_mapping] = reduce_surface(bst_surf, dipole_idx, ...
    max_distance);

% Check visually
figure;
hold on;
colors = ones(1, size(bst_surf.Vertices, 1));
cortex_patch = trimesh(bst_surf.Faces, bst_surf.Vertices(:, 1), ... 
    bst_surf.Vertices(:, 2), bst_surf.Vertices(:, 3), colors);

colors_reduced = 2*ones(1, size(bst_surf_reduced.Vertices, 1));
reduced_patch = trimesh(bst_surf_reduced.Faces, ... 
    bst_surf_reduced.Vertices(:, 1), bst_surf_reduced.Vertices(:, 2), ...
    bst_surf_reduced.Vertices(:, 3), colors_reduced);

dipole_locations_mat = bst_surf.Vertices(dipole_idx, :);
scatter3(...
    dipole_locations_mat(:, 1), ...
    dipole_locations_mat(:, 2), ...
    dipole_locations_mat(:, 3), ...
    100, 3 + (1:dipole_cnt), 'filled');
hold off;
clearvars -except bst_surf bst_surf_reduced dipole_cnt dipole_idx dipoles ...
    index_mapping max_distance;

% Flatten and visualize again
% First, flatten neigbours of the central vertex
[neighbours_idx, flat_coordinates] = flatten_around_a_vertex( ...
    bst_surf_reduced, index_mapping(dipole_idx(1)));
% Now, unfold using LSCM with above points as fixed
flat = lscm(bst_surf_reduced.Vertices, bst_surf_reduced.Faces, ...
    neighbours_idx', flat_coordinates);

% Visualize the flattened surface
figure;
hold on;
colors = ones(1, size(bst_surf_reduced.Vertices, 1));
flat_patch = trimesh(bst_surf_reduced.Faces, flat(:, 1), flat(:, 2), ...
    zeros(size(flat, 1), 1), colors);

dipole_flat_loc_mat = flat(index_mapping(dipole_idx), :);
scatter3(...
    dipole_flat_loc_mat(:, 1), ...
    dipole_flat_loc_mat(:, 2), ...
    zeros(dipole_cnt, 1), ...
    100, 3 + (1:dipole_cnt), 'filled');
hold off;
clearvars -except dipole_cnt dipole_idx index_mapping max_distance flat ...
    bst_surf_reduced;

%% Waves between peaks
% Generate waves
T = 0.12;
time_step = 0.004;
time_points = 0:time_step:T;
time_point_cnt = length(time_points);
half_width = max_distance/4;
speed = (max_distance - 2*half_width)/T;

directions = (0:7).*(pi/4);
directions(end+1) = NaN;
angular_half_width = pi/4;

basis_waves = {};
for direction = directions
    basis_waves{end+1} = generate_basis_wave( ...
            flat, bst_surf_reduced.Faces, ...
            time_points, speed, half_width, ...
            direction, angular_half_width);
end

% Visualize basis waves
figure;
hold off;
zmax = max(max([basis_waves{:}]));
for wave_id = 1:length(basis_waves)
    for t = 1:time_point_cnt
        trimesh(bst_surf_reduced.Faces, flat(:, 1), flat(:, 2), ...
            basis_waves{wave_id}(:, t));
        zlim([0, zmax]);
        pause(time_step);
    end
end

clearvars -except dipole_cnt dipole_idx index_mapping max_distance flat ...
    basis_waves bst_surf_reduced time_step;

% In bst compute sources with mne and constrained sources
% Import to Matlab as sources
% Then substitute  sources.ImageGridAmp for the basis waves
dipole_total_cnt = size(sources.ImageGridAmp, 1);
time_point_total_cnt = size(cell2mat(basis_waves), 2);
sources.ImageGridAmp = zeros(dipole_total_cnt, time_point_total_cnt);
wave_dipole_idx = find(~isnan(index_mapping));
sources.ImageGridAmp(wave_dipole_idx, :) = cell2mat(basis_waves);
sources.Time =  time_step*(1:time_point_total_cnt);
sources.Comment = 'Basis waves';
sources.DataFile = '';

clearvars -except dipole_cnt dipole_idx index_mapping max_distance flat ...
    basis_waves bst_surf_reduced sources time_step;

% Export to bst, project, import as basis_waves_projected

% Find at which time point on average the second peak is
F = basis_waves_projected.F(32:99, :);
wave_cnt = length(basis_waves);
time_point_cnt = size(F, 2)/ wave_cnt;
argmax = zeros(wave_cnt, 1);
for i = 1:wave_cnt
    time_from =  time_point_cnt * (i-1) + 1;
    time_till = time_point_cnt * i;
    [~, argmax(i)] = max(mean(abs(F(:, time_from:time_till))));
end

% It should be at time point 6, so we need to speed up our wave
speed_factor = fix(mean(argmax./5));

clearvars -except dipole_cnt dipole_idx index_mapping max_distance flat ...
    basis_waves bst_surf_reduced sources time_step speed_factor;


% Redo the waves but now speed them up and cut them at 7 time points
T = 0.12;
time_step = 1/SAMPLING_RATE;
time_points = 0:time_step:T;
time_point_cnt = length(time_points);
half_width = max_distance/4;
speed = (max_distance - 2*half_width)/T*speed_factor;

directions = (0:7).*(pi/4);
directions(end+1) = NaN;
angular_half_width = pi/4;

basis_waves = {};
for direction = directions
    basis_waves{end+1} = generate_basis_wave( ...
            flat, bst_surf_reduced.Faces, ...
            time_points, speed, half_width, ...
            direction, angular_half_width);
    basis_waves{end} = basis_waves{end}(:,1:7);
end


% Let's move new basis waves to bst
dipole_total_cnt = size(sources.ImageGridAmp, 1);
time_point_total_cnt = size(cell2mat(basis_waves), 2);
sources.ImageGridAmp = zeros(dipole_total_cnt, time_point_total_cnt);
wave_dipole_idx = find(~isnan(index_mapping));
sources.ImageGridAmp(wave_dipole_idx, :) = cell2mat(basis_waves);
sources.Time =  time_step*(1:time_point_total_cnt);
sources.Comment = 'Basis waves adjusted';
sources.DataFile = '';

% Project them on sensors in bst
% Import as basis_waves_projected
X = basis_waves_projected.F(32:99, :);
basis_wave_cnt = length(basis_waves);
% Cut into separate waves and convert each sensor-space wave matrix into a
% column
X = arrayfun(@(bw_id) X(:, (bw_id-1)*7+1:bw_id*7), 1:basis_wave_cnt, 'un', 0);
X = cell2mat(cellfun(@(wave_matrix) wave_matrix(:), X, 'un', 0));

% Import averaged spike (averaged around first peak) from bst as averaged_spike
% Unravel first 7 points of the averaged spike starting with zero
time_zero_id = find(abs(averaged_spike.Time) < 0.004/4);
Y = averaged_spike.F(32:99, time_zero_id:time_zero_id+6);
Y = Y(:);


% Do the LASSO
[B, FitInfo] = lasso(X, Y);
lassoPlot(B,FitInfo,'PlotType','Lambda','XScale','log');

% Take the first solution that has exactly one non-zero coefficient out of
% the first eight
solution_column = find(sum(B(1:8, :) ~= 0) == 1, 1);
solution_coeff = B(:, solution_column);
solution = arrayfun(@(wave, coeff) wave{:} * coeff, ...
    basis_waves, solution_coeff', 'un', 0);
solution = sum(cat(3, solution{:}), 3);

% Put the solution into sources
dipole_total_cnt = size(sources.ImageGridAmp, 1);
sources.ImageGridAmp = zeros(dipole_total_cnt, size(solution, 2));
wave_dipole_idx = find(~isnan(index_mapping));
sources.ImageGridAmp(wave_dipole_idx, :) = solution;
sources.Time =  time_step*(1:size(solution, 2));
sources.Comment = 'Solution wave';
sources.DataFile = '';

% save('161224_session.mat');

%% Find the wave to approximate activity after the second peak

% Find the point that is closest to the wave peak at time point 6
solution_directed = basis_waves{find(solution_coeff, 1)};
[~, second_peak_center_id] = max(solution_directed(:, 6));

% Visualize
figure;
hold on;
colors = ones(1, size(bst_surf_reduced.Vertices, 1));
flat_patch = trimesh(bst_surf_reduced.Faces, flat(:, 1), flat(:, 2), ...
    zeros(size(flat, 1), 1), colors);

dipole_flat_loc_mat = flat(index_mapping(dipole_idx), :);
scatter3(...
    dipole_flat_loc_mat(:, 1), ...
    dipole_flat_loc_mat(:, 2), ...
    zeros(dipole_cnt, 1), ...
    100, 3 + (1:dipole_cnt), 'filled');
scatter3(flat(second_peak_center_id, 1), ...
    flat(second_peak_center_id, 2), 0, 100, 7, 'filled');
hold off;


% Now, do the whole thing again

second_peak_center_orig_id = find(index_mapping == second_peak_center_id);
[bst_surf_reduced2, index_mapping2] = reduce_surface(bst_surf, [second_peak_center_orig_id], ...
    max_distance * 2.2);
[neighbours_idx2, flat_coordinates2] = flatten_around_a_vertex( ...
    bst_surf_reduced2, index_mapping2(second_peak_center_orig_id));

% Check visually
figure;
hold on;
colors = ones(1, size(bst_surf.Vertices, 1));
cortex_patch = trimesh(bst_surf.Faces, bst_surf.Vertices(:, 1), ... 
    bst_surf.Vertices(:, 2), bst_surf.Vertices(:, 3), colors);

colors_reduced = 2*ones(1, size(bst_surf_reduced2.Vertices, 1));
reduced_patch = trimesh(bst_surf_reduced2.Faces, ... 
    bst_surf_reduced2.Vertices(:, 1), bst_surf_reduced2.Vertices(:, 2), ...
    bst_surf_reduced2.Vertices(:, 3), colors_reduced);

hold off;


% Now, unfold using LSCM with above points as fixed
flat2 = lscm(bst_surf_reduced2.Vertices, bst_surf_reduced2.Faces, ...
    neighbours_idx2', flat_coordinates2);

% Visualize the flattened surface
figure;
hold on;
colors = ones(1, size(bst_surf_reduced2.Vertices, 1));
flat_patch = trimesh(bst_surf_reduced2.Faces, flat2(:, 1), flat2(:, 2), ...
    zeros(size(flat2, 1), 1), colors);

dipole_flat_loc_mat = flat2(index_mapping2(dipole_idx), :);
scatter3(...
    dipole_flat_loc_mat(:, 1), ...
    dipole_flat_loc_mat(:, 2), ...
    zeros(dipole_cnt, 1), ...
    100, 3 + (1:dipole_cnt), 'filled');
hold off;

basis_waves2 = {};
for direction = directions
    basis_waves2{end+1} = generate_basis_wave( ...
            flat2, bst_surf_reduced2.Faces, ...
            T*1000, speed, half_width, ...
            direction, angular_half_width, SAMPLING_RATE);
    basis_waves2{end} = basis_waves2{end}(:,1:20);
end



% Let's move new basis waves to bst
dipole_total_cnt = size(sources.ImageGridAmp, 1);
time_point_total_cnt = size(cell2mat(basis_waves2), 2);
sources.ImageGridAmp = zeros(dipole_total_cnt, time_point_total_cnt);
wave_dipole_idx2 = find(~isnan(index_mapping2));
sources.ImageGridAmp(wave_dipole_idx2, :) = cell2mat(basis_waves2);
sources.Time =  time_step*(1:time_point_total_cnt);
sources.Comment = 'Basis waves adjusted 2';
sources.DataFile = '';

% Project them on sensors in bst
% Import as basis_waves_projected
X2 = basis_waves_projected2.F(32:99, :);
basis_wave_cnt2 = length(basis_waves2);
time_point_cnt2 = size(X2, 2)/ basis_wave_cnt2;
% Cut into separate waves and convert each sensor-space wave matrix into a
% column
X2 = arrayfun(@(bw_id) X2(:, (bw_id-1)*time_point_cnt2+[1:time_point_cnt2]), ...
    1:basis_wave_cnt2, 'un', 0);
X2 = cell2mat(cellfun(@(wave_matrix) wave_matrix(:), X2, 'un', 0));

% Import averaged spike (averaged around first peak) from bst as averaged_spike
% Unravel first 7 points of the averaged spike starting with zero
time_zero_id2 = find(abs(averaged_spike2.Time) < 0.004/4);
Y2 = averaged_spike2.F(32:99, time_zero_id2 + [0:time_point_cnt2-1]);
Y2 = Y2(:);


% Do the LASSO
[B2, FitInfo2] = lasso(X2, Y2);
lassoPlot(B2,FitInfo2,'PlotType','Lambda','XScale','log');

% Take the first solution that has exactly one non-zero coefficient out of
% the first eight
solution_column2 = find(sum(B2(1:8, :) ~= 0) == 0, 1) - 1;
solution_coeff2 = B2(:, solution_column2);
solution_coeff2(end) = 0;
solution2 = arrayfun(@(wave, coeff) wave{:} * coeff, ...
    basis_waves2, solution_coeff2', 'un', 0);
solution2 = sum(cat(3, solution2{:}), 3);

% Put the solution into sources
dipole_total_cnt = size(sources.ImageGridAmp, 1);
sources.ImageGridAmp = zeros(dipole_total_cnt, size(solution2, 2));
wave_dipole_idx2 = find(~isnan(index_mapping2));
sources.ImageGridAmp(wave_dipole_idx2, :) = solution2;
sources.Time =  time_step*(1:size(solution2, 2));
sources.Comment = 'Solution wave2';
sources.DataFile = '';

save('161230.mat');

