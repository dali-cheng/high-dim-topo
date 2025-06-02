clearvars; close all; clc;

figure;
[xx, yy, zz] = sphere(); r = 0.05;
surf(r*xx, r*yy, 0.5 + r*zz, ...
    'FaceColor', 'r', 'EdgeColor', 'none'); hold on;
surf(r*xx, r*yy, -0.5 + r*zz, ...
    'FaceColor', 'b', 'EdgeColor', 'none'); hold on;

vertices = [-1, -1, -1; 1, -1, -1; 1, 1, -1; -1, 1, -1; ...
    -1, -1, 1;  1, -1, 1;  1, 1, 1;  -1, 1, 1];
faces = [1, 2, 3, 4; 5, 6, 7, 8; 1, 2, 6, 5; ...
    2, 3, 7, 6; 3, 4, 8, 7; 4, 1, 5, 8];
patch('Vertices', vertices, 'Faces', faces, 'FaceColor', [1, 1, 1] * 0.7, ...
    'FaceAlpha', 0.1, 'EdgeColor', [0, 0, 0], 'LineWidth', 1); hold on;

mArrow3([0, 0, 0], [1.5, 0, 0], 'stemWidth', 0.015, 'tipWidth', 0.05, 'FaceAlpha', 1);
mArrow3([0, 0, 0], [0, 1.5, 0], 'stemWidth', 0.015, 'tipWidth', 0.05, 'FaceAlpha', 1);
mArrow3([0, 0, 0], [0, 0, 1.5], 'stemWidth', 0.015, 'tipWidth', 0.05, 'FaceAlpha', 1);

axis equal; axis([-2, 1, -1, 1, -1, 1] * 1.5);
xlabel(''); ylabel(''); zlabel('');
grid off; axis off;
set(gca, 'fontname', 'Arial', 'fontsize', 18, 'fontweight', 'normal', ...
    'labelfontsizemultiplier', 1, 'linewidth', 2, 'Layer', 'Top', 'box', 'on');
set(gcf, 'unit', 'normalized', 'Position', [0.05 0.05 0.8 0.8], 'Color', 'w');
lightangle(138, 17); lighting gouraud;
view([117, 10.5]);