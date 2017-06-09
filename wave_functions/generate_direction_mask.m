function mask = generate_direction_mask(i, j, phi, direction, angular_half_width)
    if isnan(direction)
        mask = 1;
        return;
    end
    phi_deltas = mod((phi - direction + pi), 2*pi) - pi;
    phi_deltas(i == j) = 0;
    mask = (1 + cos(pi * phi_deltas / angular_half_width)) ...
        .* (abs(phi_deltas) <= angular_half_width);
end