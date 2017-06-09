function [W, FitInfo] = pick_waves_using_lasso(waves, meas, G, skipped_idx)
% @param waves - cell array of length L, each cell contains a wave
% a wave - NxT matrix, where
%   N - number of cortical vertices
%   T - number of time points
% @param meas - MxT matrix of measurements, where
%   M - number of sensors
% @param G - MxN matrix of the forward model
% @param skipped_idx - non-zero coefficients for these waves are not
% accounted for
%
% @returns W - vector of coefficients from LASSO
%
% Finds minimal linear combination of waves that models eeg the best
% using LASSO on unraveled projected waves as explanatory variables and
% eeg as the independent variable. Minimal - in the sense of the number of
% non-zero coefficients


    % Project waves on sensors, unravel each matrix so that each wave is
    % just one long vector, then stack those vectors side by side in an
    % (N*T)xL matrix
    unravel = @(matrix) reshape(matrix, numel(matrix), 1);
    X = cellfun(@(wave) unravel(G*wave), waves, 'un', 0);
    X = cell2mat(X);
    Y = unravel(meas);

    % Do the LASSO. Each column of B contains coefficients for one value of
    % lambda. Columns closer to the end correspond to a higher value of
    % lambda
    [B, FitInfo] = lasso(X, Y);
    FitInfo.B = B;
    
    % Return the coefficient column which is the first with the most zeros
    % not counting the waves with skipped_idx indices.
    B_ = B; B_(skipped_idx, :) = [];
    coeff_cnt = sum(abs(B_) > eps); % number of non-zero coefficients
    solution_coeff_cnt = min(nonzeros(coeff_cnt));
    solution_id = find(coeff_cnt == solution_coeff_cnt, 1);
    W = B(:, solution_id);

end