clearvars; clc; close all;
varphi_2 = -0.58 * pi;

figure;
[Omegat_plot, delta_omega_plot, spectrum1] = get_spectrum([1; 0], varphi_2);
[~, ~, spectrum2] = get_spectrum([0; 1], varphi_2);
spectrum = spectrum1 + spectrum2;
spectrum = spectrum / max(max(spectrum));

surf(Omegat_plot, delta_omega_plot, transpose(spectrum), 'EdgeColor', 'none');
colormap hot; view([0, 0, 1]); grid off;
xlim([0.4, 0.5]); ylim([-1, 1] * 0.02); yticks([-1, 0, 1] * 0.02);
xticklabels([]); yticklabels([]); 
set(gca, 'fontname', 'Arial', 'fontsize', 22, 'fontweight', 'normal', ...
    'labelfontsizemultiplier', 1, 'linewidth', 2, 'Layer', 'Top', 'Box', 'on');
set(gcf, 'unit', 'normalized', 'Position', [0.2 0.05 0.45 0.7]);


function [Omegat_plot, delta_omega_plot, spectrum] = get_spectrum(psi_in, varphi_2)
M = 3;
kappa_TR = 0.1; gamma0 = 0.02;
delta_omega_range = 0.02;

num_pts_within_roundtrip = 4000; num_roundtrips = 201;
spectrum = zeros(num_pts_within_roundtrip * num_roundtrips, 1);
Omegat_list = linspace(0, num_roundtrips, numel(spectrum) + 1) * 2*pi; Omegat_list(end) = [];
delta_omega_list = linspace(-1, 1, numel(spectrum) + 1) * delta_omega_range * 2*pi; delta_omega_list(end) = [];
[Omegat_plot, delta_omega_plot] = meshgrid(linspace(-0.5, 0.5, num_pts_within_roundtrip), ...
    linspace(-1, 1, num_roundtrips - 1) * delta_omega_range);

T = @(x) transmittance(x, M, kappa_TR, varphi_2);

for loop_index = 1 : numel(spectrum)
    Omegat = Omegat_list(loop_index);
    delta_omega = delta_omega_list(loop_index);
    P = exp(1i * delta_omega) * exp(-gamma0);
    psi_ss = inv(eye(2) - T(Omegat) * P) * psi_in;
    spectrum(loop_index) = psi_ss' * psi_ss;
end
spectrum(1 : num_pts_within_roundtrip/2) = [];
spectrum(end - num_pts_within_roundtrip/2 + 1 : end) = [];
spectrum = reshape(spectrum, num_pts_within_roundtrip, num_roundtrips - 1);
end


function y = transmittance(Omegat, M, kappa_TR, varphi_2)
U1 = [-1, 1; 1, 1] / sqrt(2);
U2 = [1i, -1i; 1, 1] / sqrt(2);

k1 = Omegat; k2 = M * Omegat + varphi_2;

y1 = [exp(-1i * kappa_TR * (3 * cos(k1) + 4 * cos(k2))) * exp(-0.2 * kappa_TR * cos(k1)), 0;
      0, exp(+1i * kappa_TR * (3 * cos(k1) + 4 * cos(k2))) * exp(+0.2 * kappa_TR * cos(k1))];
y2 = [exp(-1i * kappa_TR * sin(k1)) * exp(-0.2 * kappa_TR * sin(k1)), 0;
      0, exp(+1i * kappa_TR * sin(k1)) * exp(+0.2 * kappa_TR * sin(k1))];

y = (U1 * y1 * U1' + U2 * y2 * U2')/2;
end