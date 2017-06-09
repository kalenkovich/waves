test_mesh.Vertices = [
    0,  0;
   -1,  0;
    0,  1;   
    1,  0;
    0, -1;
   -1,  1;
    1,  1;
    1, -1;
   -1, -1
    ];
test_mesh.Vertices(:, 3) = 0;

test_mesh.Faces = [
    1, 2, 3;
    1, 3, 4;
    1, 4, 5;
    1, 5, 2;
    2, 3, 6;
    3, 4, 7;
    4, 5, 8;
    5, 2, 9
    ];

PARAMS.max_distance = 5;
test_mesh = calculate_distances_and_angles(test_mesh, PARAMS);

labels = cellstr(num2str([1:9]'));
dx = 0.1; dy = 0.1; % displacement so the text does not overlay the data points


figure;
hold on;
x = test_mesh.Vertices(:, 1);
y = test_mesh.Vertices(:, 2);
trimesh(test_mesh.Faces, x, y);
scatter(x, y, 10, 1:size(test_mesh.Vertices, 1), 'filled');
text(x + dx, y + dy, labels);

figure;
[x, y] = pol2cart(test_mesh.angles(1,:), test_mesh.distances(1, :));
hold on;
trimesh(test_mesh.Faces, x, y);
scatter(x, y, 10, 1:size(test_mesh.Vertices, 1), 'filled');
text(x + dx, y + dy, labels);