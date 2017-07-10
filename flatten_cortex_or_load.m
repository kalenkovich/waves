function cortex = flatten_cortex_or_load(PARAMS, resolution)
    processed_cortex_fname = sprintf('%s_%s_15K_%dcm_processed.mat', ...
        PARAMS.cortex.(resolution).comment, PARAMS.subject_name, ...
        PARAMS.max_distance*100);
    cortex_comment = PARAMS.cortex.low.comment;

    if ~exist(processed_cortex_fname, 'file')
        cortex = from_bst_get_surface(cortex_comment, PARAMS);
        disp('Flattening the cortex. This will take several hours.');
        cortex = calculate_distances_and_angles(cortex, PARAMS);
        save(processed_cortex_fname, '-struct', 'cortex');
    else
        disp('Loading precalculated flattened cortex');
        cortex = load(processed_cortex_fname);
    end
end