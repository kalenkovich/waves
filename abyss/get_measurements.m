function meas = get_measurements(bst_avg, time_period_in_ms, SAMPLING_RATE, ...
                 MEG_CHANNELS_IDX)
    SPIKE2_COMMENT = 'Avg: spike_peak2 (8)';
    if strcmpi(bst_avg.Comment, SPIKE2_COMMENT)
        sample_cnt = ceil(SAMPLING_RATE/1000*time_period_in_ms);
        time_zero_point = find(abs(bst_avg.Time) < 1/SAMPLING_RATE/10);
        meas = bst_avg.F(MEG_CHANNELS_IDX, ...
            time_zero_point + [0:sample_cnt]);
    else
        error(['Smth wrong, probaly ' SPIKE2_COMMENT ...
               ' does not exist']);
    end
end