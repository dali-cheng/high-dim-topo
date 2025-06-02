clearvars; close all; clc;

num_k1 = 61; num_k3 = 61;
k1_list = linspace(-pi, pi, num_k1);
k3_list = linspace(-pi, pi, num_k3);
E = zeros(num_k1, num_k3, 2);

for k1_index = 1 : num_k1
    for k3_index = 1 : num_k3
        k1 = k1_list(k1_index); k3 = k3_list(k3_index);
        H = sin(k1) * pauli(1) + (1 - cos(k1) - cos(k3)) * pauli(3);
        E(k1_index, k3_index, :) = sort(eig(H));
    end
end

[k1_plot, k3_plot] = meshgrid(k1_list, k3_list);
figure;
mesh(k1_plot/pi, k3_plot/pi, transpose(squeeze(E(:, :, 1))), ...
     'EdgeColor', 'k'); hold on;
mesh(k1_plot/pi, k3_plot/pi, transpose(squeeze(E(:, :, 2))), ...
     'EdgeColor', 'k'); hold off;

xlabel(''); xticklabels(''); xticks([-1, 0, 1]); 
ylabel(''); yticklabels(''); yticks([-1, 0, 1]); 
zlabel(''); zticklabels(''); zticks([-3, 0, 3]); 
axis([-1, 1, -1, 1, -3, 3]); view([-24, 8]);
grid off;
set(gca, 'fontname', 'Arial', 'fontsize', 30, 'fontweight', 'normal', ...
    'labelfontsizemultiplier', 1, 'linewidth', 3, 'Layer', 'Top', 'box', 'on');
set(gcf, 'unit', 'normalized', 'Position', [0.2 0.05 0.5 0.6]);
view([144, 6.5]);