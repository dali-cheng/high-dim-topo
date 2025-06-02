clearvars; close all; clc;
M = 3; L = 21;

H = kron(diag(ones(1, M^2*L-1), -1),     1i * pauli(1) - pauli(3)) + ...
    kron(diag(ones(1, M^2*L-M), -M),     -pauli(3)) + ...
    kron(diag(ones(1, M^2*L-M^2), -M^2), 1i * pauli(2) - pauli(3));
H = (H + H' + kron(diag(ones(1, M^2*L)), 4 * pauli(3))) / 2;

delta_omega = 0; eta_0 = 0.01;
a_in = zeros(2*M^2*L, 1); a_in(1) = 1;
a_ss = inv((delta_omega - 0.5i * eta_0) * eye(size(H)) - H) * a_in;
site_population = zeros(M^2*L, 1);
for l = 1 : M^2*L
    site_population(l) = sum(abs(a_ss(2*l-1 : 2*l)).^2);
end


site_radius = site_population / max(site_population) * 0.48 + 0.02;
% Convert to [0.02, 0.5] for visualization purposes
[xx, yy, zz] = sphere(50);
figure; my_color = hot;
for l = 1 : M^2*5
    [x, y, z] = get_site_index(l, M);
    surf(x + site_radius(l) * xx, ...
         y + site_radius(l) * yy, ...
         z + site_radius(l) * zz, ...
         'FaceColor', my_color(round((site_radius(l) - 0.02) / 0.48 * 240 + 1), :), ...
         'EdgeColor', 'none'); hold on;
end
grid off; axis equal; view([-130, 12]);
lightangle(220, 20); lighting gouraud; material dull;
xlim([-0.5, M - 0.5]); ylim([-0.5, M - 0.5]); zlim([-0.5, 5 - 0.5]); 
xticks([0, M-1]); yticks([0, M-1]); zticks([0, 5-1]);
xticklabels([]); yticklabels([]); zticklabels([]);
set(gca, 'fontname', 'Arial', 'fontsize', 22, 'fontweight', 'normal', ...
    'labelfontsizemultiplier', 1, 'linewidth', 2, 'Layer', 'Top', 'Box', 'on');
set(gcf, 'unit', 'normalized', 'Position', [0.2 0.05 0.4 0.8]);


layer_population = zeros(L, 1);
for l = 1 : L
    layer_population(l) = sum(site_population(M^2*(l-1)+1 : M^2*l));
end
layer_population = layer_population / sum(layer_population);

figure; 
semilogy(0:L-1, layer_population, 'LineWidth', 2, 'Color', 'k');
view(90, -90);
xlim([0, L-1]); ylim([1e-10, 1]); xticks([0, L-1]);
xticklabels([]); yticklabels([]);
set(gca, 'fontname', 'Arial', 'fontsize', 22, 'fontweight', 'normal', ...
    'labelfontsizemultiplier', 1, 'linewidth', 2, 'Layer', 'Top', 'Box', 'on');
set(gcf, 'unit', 'normalized', 'Position', [0.2 0.05 0.25 0.8]);


function [x, y, z] = get_site_index(l, M)
l = l - 1;
z = floor(l / M^2);
y = floor((l - z * M^2) / M);
x = l - z * M^2 - y * M;
end