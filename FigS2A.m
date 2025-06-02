clearvars; close all; clc;

k1_list = linspace(0.8, 1, 1001) * pi;
k2_list = linspace(0.1, 0.3, 1001) * pi;
[k1_plot, k2_plot] = meshgrid(k1_list, k2_list);
real_E = zeros(numel(k1_list), numel(k2_list), 2);
imag_E = zeros(numel(k1_list), numel(k2_list), 2);

for k1_index = 1 : numel(k1_list)
    for k2_index = 1 : numel(k2_list)
        k1 = k1_list(k1_index);
        k2 = k2_list(k2_index);

        H = ((3 - 0.2i) * cos(k1) + 4 * cos(k2)) * pauli(1) ...
            + (1 - 0.2i) * sin(k1) * pauli(2);

        eig_val = eig(H);
        real_E(k1_index, k2_index, :) = real(eig_val);
        imag_E(k1_index, k2_index, :) = imag(eig_val);
    end
end

figure;
mesh(k1_plot/pi, k2_plot/pi, transpose(squeeze(real_E(:, :, 1))), ...
    'FaceColor', 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.15, 'LineWidth', 1); hold on;
mesh(k1_plot/pi, k2_plot/pi, transpose(squeeze(real_E(:, :, 2))), ...
    'FaceColor', 'b', 'EdgeColor', 'none', 'FaceAlpha', 0.15, 'LineWidth', 1);
xlim([0.8, 1]); xticks([0.8, 1]); xticklabels([]);
ylim([0.1, 0.3]); yticks([0.1, 0.3]); yticklabels([]);
zlim([-1, 1] * 2); zticks([-1, 1] * 2); zticklabels([]);
view([-146, 6]);
set(gca, 'fontname', 'Arial', 'fontsize', 22, 'fontweight', 'normal', ...
    'labelfontsizemultiplier', 1, 'linewidth', 2, 'Layer', 'Top', 'Box', 'on');
set(gcf, 'unit', 'normalized', 'Position', [0.2 0.05 0.45 0.7]);