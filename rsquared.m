v_id = randi(size(sensor_waves, 1));
s_id = randi(size(sensor_waves, 2));
d_id = randi(size(sensor_waves, 3));

W = squeeze(sensor_waves(v_id, s_id, d_id, :, :));

Y = meas{1}.F(:, 1:6);
w = repmat(sqrt(sum(Y.^2, 1)), 68, 1);

x = lscov(W(:), Y(:), 1./w(:).^2);

norm(x*W(:)./w(:))/norm(Y(:)./w(:))


scores(v_id, s_id, d_id, 1)