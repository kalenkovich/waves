function ndim = cell2ndim(cell_arr)
% Converts M-dimensional cell array of N-dimensional numerical arrays into
% (M+N)-dimensional numerical array. M is inferred from size(cell_arr) and
% so can never be less than 1.
    M = length(size(cell_arr));
    M_ones = num2cell(ones(M, 1));
    num_arr_sz = num2cell(size(cell_arr{1}));
    cell_arr = cellfun(@(num_arr) ...
        reshape(num_arr, M_ones{:}, num_arr_sz{:}), cell_arr, 'un', 0);
    ndim = cell2mat(cell_arr);
end