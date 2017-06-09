function gain_matrix = from_bst_get_gain_matrix(comment, PARAMS)

    % Initialize Brainstorm and load the right protocol
    from_bst_initialize(PARAMS);
    
    % Find the head model
    subject_name = PARAMS.subject_name;
    study_name = PARAMS.study_name;
    subject_struct = bst_get('Subject', subject_name);
    studies_with_subject = bst_get('StudyWithSubject', subject_struct.FileName);
    [~, iStudy] = ismember(study_name, {studies_with_subject.Name});
    study_info = studies_with_subject(iStudy);
    
    head_model_id = find(...
        cellfun(@(Comment) strcmpi(Comment, comment), ...
            {study_info.HeadModel.Comment}));
        
    % Check that exactly one was found
    assert(length(head_model_id) == 1, ...
        '%d head models matched. Not good', length(head_model_id));    
        
    % Load the head model
    head_model_relative_path = study_info.HeadModel(head_model_id).FileName;
    bst_head_model = load(file_fullpath(head_model_relative_path));
    
    % Convert the gain matrix to fixed orientations of dipoles. Then take
    % only the rows corresponding to the channels we need.
    channels_idx = PARAMS.channels_idx;
    gain_matrix = bst_gain_orient(bst_head_model.Gain, ...
        bst_head_model.GridOrient);
    gain_matrix = gain_matrix(channels_idx, :);

end