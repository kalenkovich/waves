G = from_bst_get_gain_matrix(PARAMS.forward.name, PARAMS);
T = get_time_point_cnt(PARAMS);
Tw = 2*T;
N = 500;
Fs = PARAMS.sampling_rate;

sim_cnt = 1000;
noise = arrayfun(@(i) GenerateBrainNoise(G, T, Tw, N, Fs), 1:sim_cnt, 'un', 0);
noise = cellfun(@(x) reshape(x, 1, size(x, 1), size(x, 2)), noise, 'un', 0);
noise = cell2mat(noise');
% (simulation_number, sensor_number, time_point) 3D array

noise_power = sum(sum(noise.^2, 2), 3);
noise = bsxfun(@rdivide, noise, sqrt(noise_power));

% waves from one vertex at some speeds
v0 = randi(size(sensor_waves, 1));
speed_ids = [1, 3, 9];
speed_cnt = length(speed_ids);
direction_cnt = size(sensor_waves, 3);
sensor_cnt = size(sensor_waves, 4);
waves = squeeze(sensor_waves(v0, speed_ids, :, :, :));
% (speed, direction, sensor, time)

% One wave that will represent the signal
s0 = randi(3);
d0 = randi(8) + 1; % 0th direction corresponds to no direction
the_wave = squeeze(waves(s0, d0, :, :));
the_wave_power = sum(sum(the_wave.^2));
the_wave = the_wave ./ sqrt(the_wave_power);
the_wave = reshape(the_wave, 1, size(the_wave, 1), size(the_wave, 2));

% SNR = summed squared magnitude of the signal due to wave over summed
% squared magnitude of noise
SNRs = [-7:-5];
results = cell(size(SNRs));
for i = 1:length(SNRs)
    SNR = SNRs(i);
    power_ratio = 10^(SNR/10);
    the_wave_with_noise = bsxfun(@plus, the_wave * sqrt(power_ratio), noise);
    
    Y = the_wave_with_noise(2, :)';
    X = reshape(waves, speed_cnt*direction_cnt, ...
                sensor_cnt * T)';
    [B,FitInfo] = lasso(X,Y);
    lambdas = FitInfo.Lambda;
    lambdas = linspace(lambdas(95), lambdas(end));
    lambda_cnt = length(lambdas);
    B_all = NaN(sim_cnt, size(B, 1), size(B, 2));
    for j = 1:sim_cnt
        Y = the_wave_with_noise(j, :)';
        B_all(j, :, :) = lasso(X, Y);
    end
    B_all = reshape(B_all, sim_cnt, speed_cnt, direction_cnt, lambda_cnt);
    results{i} = B_all;
end

got_the_right_wave = cellfun(...
    @(B_all) squeeze(B_all(:, s0, d0, :)~=0), results, 'un', 0);


directed_wave_cnt = cellfun(...
    @(B_all) squeeze(sum(...
        reshape(B_all(:, :, 2:end, :)~=0, sim_cnt, [], lambda_cnt), 2)), ...
    results, 'un', 0);

got_a_wrong_wave = cellfun(...
    @(a, b) a > b, directed_wave_cnt, got_the_right_wave, 'un', 0);

wrong_wave_cnt = cellfun(...
    @(a, b) a - b, directed_wave_cnt, got_the_right_wave, 'un', 0);

sensitivity = cellfun(@(x) sum(x, 1)'/sim_cnt, got_the_right_wave, 'un', 0);
specificity = cellfun(@(x) 1- sum(x, 1)'/sim_cnt, got_a_wrong_wave, 'un', 0);
specificity2 = cellfun(@(x) ...
    1- sum(x, 1)'/sim_cnt/((direction_cnt-1)*speed_cnt - 1), ...
    wrong_wave_cnt, 'un', 0);

[sensitivity{1} specificity{1} specificity2{1}]

figure;
for i = [1:3]
    hold on;
    subplot(1, 3, i);
    line(1-specificity2{i}, sensitivity{i});
    line([0, 1], [0, 1]);
    hold off;
end

figure;
for i = [1:3]
    hold on;
    subplot(1, 3, i);
    
    line(1-specificity2{i}, sensitivity{i});
    line([0, 1], [0, 1]);
    
    axis square;
    set(gca, 'xtick',[0 1]);
    set(gca, 'ytick',[0 1]);
    title(sprintf('SNR = %d', SNRs(i)));
    xlabel('1-specificity');
    
    hold off;
end

subplot(1, 3, 1);
ylabel('sensitivity');