function for_evgenii2

load for_evgenii2

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

for ii = 1 : length(distances_LOS_FLOYD)

  % print the distance
  fprintf('LOS distance from vertex %3d to vertex %3d = %.3f\n', ...
          vertex_pairs(ii,1), vertex_pairs(ii,2), distances_LOS(ii));

  % plot the path
  if ~isempty(paths_LOS{ii})
    ph_LOS = plot3(paths_LOS{ii}.X, ...
                   paths_LOS{ii}.Y, ...
                   paths_LOS{ii}.Z, 'g');
    set(ph_LOS, 'linewidth', 2);
  end
  
  % print the distance
  fprintf('LOS-FLOYD distance from vertex %3d to vertex %3d = %.3f\n', ...
          vertex_pairs(ii,1), vertex_pairs(ii,2), distances_LOS_FLOYD(ii));

  % plot the path
  ph_LOS_FLOYD = plot3(paths_LOS_FLOYD{ii}.X, ...
                       paths_LOS_FLOYD{ii}.Y, ...
                       paths_LOS_FLOYD{ii}.Z, 'r');
  set(ph_LOS_FLOYD, 'linewidth', 2);
  
  fprintf('\n');
  
  pause
  if ~isempty(paths_LOS{ii})
    delete(ph_LOS)
  end
  delete(ph_LOS_FLOYD)
  
end

return
