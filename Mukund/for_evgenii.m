function for_evgenii

load for_evgenii

% show the patch
figure(1); clf
ph = patch('Vertices', vertices, 'Faces', faces, ...
           'FaceColor', [0.6 0.6 0.6], ...
           'LineWidth', 1.0); 
hold on
axis off
axis vis3d 
axis equal
lighting gouraud
material dull
view(+150, 15);
camlight(+40,30);

for ii = 1 : length(distances)

  % print the distance
  fprintf('distance from vertex %3d to vertex %3d = %.2f\n', ...
          vertex_pairs(ii,1), vertex_pairs(ii,2), distances(ii));

  % plot the path
  ph = plot3(paths(ii).X, ...
             paths(ii).Y, ...
             paths(ii).Z, 'r');
  set(ph, 'linewidth', 2);
  pause
  delete(ph)
  
end

return
