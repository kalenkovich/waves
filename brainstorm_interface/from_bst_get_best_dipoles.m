function best_vertex_idx = from_bst_get_best_dipoles(meas, ...
    dipole_kernel, vertex_cnt, PARAMS)

    perf = abs(dipole_kernel * meas);
    [~, sortIndex] = sort(perf, 'descend'); 
    best_vertex_idx = sortIndex(1:vertex_cnt); 
    
end