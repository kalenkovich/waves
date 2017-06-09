function surface = from_bst_get_surface(comment, PARAMS)

    % Initialize Brainstorm and load the right protocol
    from_bst_initialize(PARAMS);
    
    % Find the right surface
    subject_info = bst_get('Subject', PARAMS.subject_name);
    surface_id = find(...
        cellfun(@(Comment) strcmpi(Comment, comment), ...
            {subject_info.Surface.Comment}));
        
    % Check that exactly one was found
    assert(length(surface_id) == 1, ...
        '%d head models matched. Not good', length(surface_id));     
    
    % Load the surface
    surface_relative_path = subject_info.Surface(surface_id).FileName;
    surface = load(file_fullpath(surface_relative_path));

end