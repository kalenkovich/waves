function time_points = get_time_points(PARAMS)
    sampling_rate = PARAMS.sampling_rate;
    duration = PARAMS.wave.duration;
    time_points = 0:1/sampling_rate:duration;
end