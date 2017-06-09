function to_bst_wave(wave, comment, time_points, PARAMS)

    sources_struct = load('sources_stub.mat');
    sources_struct.ImageGridAmp = wave;
    time_point_cnt = size(wave, 2);
    if isempty(time_points)
        sources_struct.Time = [0:time_point_cnt-1] / PARAMS.sampling_rate;
    else
        sources_struct.Time = time_points;
    end
    sources_struct.Comment = comment;
    
    [~, iStudy] = bst_get('Study');
    db_add(iStudy, sources_struct);
end