function dipoles = from_bst_get_dipoles(comment, PARAMS)

    % Initialize Brainstorm and load the right protocol
    from_bst_initialize(PARAMS);
    
    % Find the head model
    study_info = bst_get('Study');
    dipoles_id = find(...
        cellfun(@(Comment) strcmpi(Comment, comment), ...
            {study_info.Dipoles.Comment}));
    
    % Check that exactly one was found
    assert(length(dipoles_id) == 1, ...
        '%d dipole files matched. Not good', length(dipoles_id));    
        
    % Load the head model
    dipoles_relative_path = study_info.Dipoles(dipoles_id).FileName;
    dipoles = load(file_fullpath(dipoles_relative_path));

end