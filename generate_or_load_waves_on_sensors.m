function sensor_waves = generate_or_load_waves_on_sensors(cortex, PARAMS, resolution)
    waves_fname = sprintf('%s_%s_15K_%dcm_waves.mat', ...
        PARAMS.cortex.(resolution).comment, PARAMS.subject_name, ...
        PARAMS.max_distance*100);

    if ~exist(waves_fname, 'file')
        disp('Generating waves. This will take several minutes.');
        sensor_waves = generate_waves_on_sensors(cortex, PARAMS);
        disp('Saving waves. This will take a couple minutes.');
        save(waves_fname, 'PARAMS', 'resolution', 'sensor_waves', '-v7.3');
    else
        disp('Loading precalculated waves. This will take a minute.');
        sensor_waves = load(waves_fname);
        sensor_waves = sensor_waves.sensor_waves;
    end
end