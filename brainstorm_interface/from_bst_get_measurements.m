function meas = from_bst_get_measurements(comment, PARAMS)
% Find measurements in Brainstorm whose name starts with comment

    % Initialize Brainstorm and load the right protocol
    from_bst_initialize(PARAMS);
    
    % Find the right file(s)
    file_info = bst_process('CallProcess', 'process_select_files_data', [], [], ...
        'subjectname',   PARAMS.subject_name, ...
        'tag',           comment);
    
    % Check that at least one was found
    assert(length(file_info) >= 1, ...
        'No files matched the %s tag/comment. Not good. Not good at all.', comment);
    
    % Brainstorm function finds all measurements whose name *contains*
    % comment while we only need those starting with it
    starts_with_comment_bool = strncmpi(comment, {file_info.Comment}, length(comment));
    file_info(~starts_with_comment_bool) = [];
    
    fprintf(['Found %d files with measurements whose ''comment''' ...
        ' field started with ''%s''\n'], length(file_info), comment);
    
    % Load the files
    meas = cellfun(@(fname) load(file_fullpath(fname)), ...
        {file_info.FileName}, 'un', 0);
    
    % We only need some of the channels and we only need the data itself
    % and the time points
    for i = 1:length(meas)
        meas_tmp.F = meas{i}.F(PARAMS.channels_idx, :);
        meas_tmp.Time = meas{i}.Time;
        meas{i} = meas_tmp;
    end
    
end