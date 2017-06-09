function norms = norms_of_rows(matrix)
    norms = sqrt(sum(matrix.^2, 2));
end