clearvars; close all; clc;

k3 = linspace(-pi/2, pi/2, 400);
k5 = acos(1 - cos(k3));
k3 = [k3, flip(k3)]; k5 = [k5, -flip(k5)];

figure;
plot3(zeros(size(k3)), k3/pi, k5/pi, ...
    'Color', 'r', 'LineWidth', 3); hold on;

b = 0.5;
k1 = linspace(-asin(b), asin(b), 400);
k5_1 = acos(1 - cos(k1) - sqrt(b^2 - (sin(k1)).^2));
k5_2 = -acos(1 - cos(k1) - sqrt(b^2 - (sin(k1)).^2));
k5_3 = acos(1 - cos(k1) + sqrt(b^2 - (sin(k1)).^2));
k5_4 = -acos(1 - cos(k1) + sqrt(b^2 - (sin(k1)).^2));

k1 = [k1, flip(k1)];
plot3(k1/pi, zeros(size(k1)), [k5_1, flip(k5_3)]/pi, ...
    'Color', 'b', 'LineWidth', 3); hold on;
plot3(k1/pi, zeros(size(k1)), [k5_2, flip(k5_4)]/pi, ...
    'Color', 'b', 'LineWidth', 3); hold on;

fill3([0, 0, 0, 0], [-1, -1, 1, 1], [-1, 1, 1, -1], ...
    [1, 0.9, 0.9], 'FaceAlpha', 0.5, 'EdgeColor', 'none'); hold on;
fill3([-1, 1, 1, -1], [0, 0, 0, 0], [-1, -1, 1, 1], ...
    [0.9, 0.9, 1], 'FaceAlpha', 0.5, 'EdgeColor', 'none');

axis equal; grid off; view(157, 9.9);
axis([-1, 1, -1, 1, -1, 1]);
xlabel(''); xticklabels([]); xticks([-1, 0, 1]);
ylabel(''); yticklabels([]); yticks([-1, 0, 1]);
zlabel(''); zticklabels([]); zticks([-1, 0, 1]);
set(gca, 'fontname', 'Arial', 'fontsize', 18, 'fontweight', 'normal', ...
    'labelfontsizemultiplier', 1, 'linewidth', 3, 'Layer', 'Top', 'box', 'on');
set(gcf, 'unit', 'normalized', 'Position', [0.05 0.05 0.8 0.8], 'Color', 'w');