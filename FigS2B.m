clearvars; close all; clc;
varphi_2 = -0.58 * pi;

k1_list = linspace(0.8, 1, 2001) * pi;
real_E = zeros(numel(k1_list), 2);
imag_E = zeros(numel(k1_list), 2);

for k1_index = 1 : numel(k1_list)
    k1 = k1_list(k1_index);
    k2 = 3 * k1 + varphi_2;

    H = ((3 - 0.2i) * cos(k1) + 4 * cos(k2)) * pauli(1) ...
            + (1 - 0.2i) * sin(k1) * pauli(2);

    eig_val = eig(H);
    real_E(k1_index, :) = sort(real(eig_val));
    imag_E(k1_index, :) = sort(imag(eig_val));
end

figure;
plot(k1_list/pi, real_E(:, 1), "Color", "k", "LineWidth", 2); hold on;
plot(k1_list/pi, real_E(:, 2), "Color", "k", "LineWidth", 2);
xlim([min(k1_list)/pi, max(k1_list)/pi]); xticks([min(k1_list)/pi, max(k1_list)/pi]); 
ylim([-1, 1] * 2); yticks([-1, 1] * 2); 
xticklabels([]); yticklabels([]);
set(gca, 'fontname', 'Arial', 'fontsize', 22, 'fontweight', 'normal', ...
    'labelfontsizemultiplier', 1, 'linewidth', 2, 'Layer', 'Top', 'Box', 'on');
set(gcf, 'unit', 'normalized', 'Position', [0.2 0.05 0.45 0.7]);