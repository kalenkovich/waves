function wave = generate_ripple_wave(r, speed, t, sigma)    
    mu = speed * t;
    wave = (1 + cos(pi * (r - mu) / sigma)) .* (abs(r - mu) <= sigma); 
end