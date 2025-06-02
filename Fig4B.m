clearvars; clc; close all;

u = 0.0; b = 0.5;
varphi_2 = 0;
varphi_3 = 0;
varphi_4 = 0;
varphi_5 = (1/3) * pi;

figure;
[Omegat_plot, delta_omega_plot, spectrum1] = get_spectrum( ...
    [1; 0; 0; 0], u, b, varphi_2, varphi_3, varphi_4, varphi_5);
[~, ~, spectrum2] = get_spectrum( ...
    [0; 1; 0; 0], u, b, varphi_2, varphi_3, varphi_4, varphi_5);
[~, ~, spectrum3] = get_spectrum( ...
    [0; 0; 1; 0], u, b, varphi_2, varphi_3, varphi_4, varphi_5);
[~, ~, spectrum4] = get_spectrum( ...
    [0; 0; 0; 1], u, b, varphi_2, varphi_3, varphi_4, varphi_5);
spectrum = spectrum1 + spectrum2 + spectrum3 + spectrum4;
spectrum = spectrum / max(max(spectrum));


surf(Omegat_plot, delta_omega_plot, transpose(spectrum), 'EdgeColor', 'none');
colormap hot; view([0, 0, 1]); grid off;
xlim([-0.025, 0.025]); xticks([-0.025, 0, 0.025]); clim([0, 1]);
ylim([-1, 1] * 0.015); yticks([-1, 0, 1] * 0.015);
xticklabels([]); yticklabels([]); pbaspect([3, 4, 1]);
set(gca, 'fontname', 'Arial', 'fontsize', 22, 'fontweight', 'normal', ...
    'labelfontsizemultiplier', 1, 'linewidth', 2, 'Layer', 'Top', 'Box', 'on');
set(gcf, 'unit', 'normalized', 'Position', [0.2 0.05 0.4 0.8]);


function [Omegat_plot, delta_omega_plot, spectrum] = get_spectrum( ...
    psi_in, u, b, varphi_2, varphi_3, varphi_4, varphi_5)
M = 3;
kappa_TR = 0.1; gamma0 = 0.003;
delta_omega_range = 0.015;

num_pts_within_roundtrip = 5000; num_roundtrips = 201;
spectrum = zeros(num_pts_within_roundtrip * num_roundtrips, 1);
Omegat_list = linspace(0, num_roundtrips, numel(spectrum) + 1) * 2*pi; Omegat_list(end) = [];
delta_omega_list = linspace(-1, 1, numel(spectrum) + 1) * delta_omega_range * 2*pi; delta_omega_list(end) = [];
[Omegat_plot, delta_omega_plot] = meshgrid(linspace(-0.5, 0.5, num_pts_within_roundtrip), ...
    linspace(-1, 1, num_roundtrips - 1) * delta_omega_range);

T = @(x) transmittance(x, M, kappa_TR, varphi_2, varphi_3, varphi_4, varphi_5, u, b);

for loop_index = 1 : numel(spectrum)
    Omegat = Omegat_list(loop_index);
    delta_omega = delta_omega_list(loop_index);
    P = exp(1i * delta_omega) * exp(-gamma0);
    psi_ss = inv(eye(4) - T(Omegat) * P) * psi_in;
    spectrum(loop_index) = psi_ss' * psi_ss;
end
spectrum(1 : num_pts_within_roundtrip/2) = [];
spectrum(end - num_pts_within_roundtrip/2 + 1 : end) = [];
spectrum = reshape(spectrum, num_pts_within_roundtrip, num_roundtrips - 1);
end

function y = transmittance(Omegat, M, kappa_TR, varphi_2, varphi_3, varphi_4, varphi_5, u, b)
U1 = [1, 0, 0, -1; 1, 0, 0, 1; 0, -1, -1, 0; 0, 1, -1, 0]/sqrt(2);
U2 = [-1i, 0, 0, -1i; 1, 0, 0, -1; 0, 1i, 1i, 0; 0, 1, -1, 0]/sqrt(2);
U3 = [0, 1, 0, 0; 0, 0, 0, 1; 0, 0, 1, 0; 1, 0, 0, 0];
U4 = [-1, 0, 0, -1; 0, 1, 1, 0; -1, 0, 0, 1; 0, 1, -1, 0]/sqrt(2);
U5 = [-1, 0, 0, -1; 0, -1i, -1i, 0; -1i, 0, 0, 1i; 0, 1, -1, 0]/sqrt(2);
U6 = [1, 0, 0, -1; 0, -1i, -1i, 0; -1i, 0, 0, -1i; 0, 1, -1, 0]/sqrt(2);

k1 = Omegat;
k2 = M * Omegat + varphi_2;
k3 = M^2 * Omegat + varphi_3;
k4 = M^3 * Omegat + varphi_4;
k5 = M^4 * Omegat + varphi_5;

y1 = get_modulation_matrix(sin(k1), kappa_TR);
y2 = get_modulation_matrix(sin(k2), kappa_TR);
y3 = get_modulation_matrix(sin(k3), kappa_TR);
y4 = get_modulation_matrix(sin(k4), kappa_TR);
y5 = get_modulation_matrix( ...
    4 + u - cos(k1) - cos(k2) - cos(k3) - cos(k4) - cos(k5), kappa_TR);
y6 = get_modulation_matrix(b, kappa_TR);

y = (U1 * y1 * U1' + U2 * y2 * U2' + U3 * y3 * U3' + ...
     U4 * y4 * U4' + U5 * y5 * U5' + U6 * y6 * U6')/6;
end

function y = get_modulation_matrix(x, kappa_TR)
y = diag([exp(+1i * kappa_TR * x), exp(+1i * kappa_TR * x), ...
          exp(-1i * kappa_TR * x), exp(-1i * kappa_TR * x)]);
end