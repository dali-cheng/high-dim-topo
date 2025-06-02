clearvars; clc; close all;
u = 0;
varphi_2 = 0.0 * pi;
varphi_3 = (1/2) * pi;

figure;
[Omegat_plot, delta_omega_plot, spectrum1] = get_spectrum([1; 0], u, varphi_2, varphi_3);
[~, ~, spectrum2] = get_spectrum([0; 1], u, varphi_2, varphi_3);
spectrum = spectrum1 + spectrum2;
spectrum = spectrum / max(max(spectrum));

%% Plotting
surf(Omegat_plot, delta_omega_plot, transpose(spectrum), 'EdgeColor', 'none');
colormap hot; view([0, 0, 1]); grid off;
xlim([-0.5, 0.5]); xticks([-0.5, 0, 0.5]); clim([0, 0.5]);
ylim([-1, 1] * 0.03); yticks([-1, 0, 1] * 0.03);

xticklabels([]); yticklabels([]); pbaspect([3, 4, 1]);
set(gca, 'fontname', 'Arial', 'fontsize', 22, 'fontweight', 'normal', ...
    'labelfontsizemultiplier', 1, 'linewidth', 2, 'Layer', 'Top', 'Box', 'on');
set(gcf, 'unit', 'normalized', 'Position', [0.2 0.05 0.4 0.8]);


function [Omegat_plot, delta_omega_plot, spectrum] = get_spectrum(psi_in, u, varphi_2, varphi_3)
M = 3;
kappa_TR = 0.1; gamma0 = 0.008;
delta_omega_range = 0.03;

num_pts_within_roundtrip = 400; num_roundtrips = 201;
spectrum = zeros(num_pts_within_roundtrip * num_roundtrips, 1);
Omegat_list = linspace(0, num_roundtrips, numel(spectrum) + 1) * 2*pi; Omegat_list(end) = [];
delta_omega_list = linspace(-1, 1, numel(spectrum) + 1) * delta_omega_range * 2*pi; delta_omega_list(end) = [];
[Omegat_plot, delta_omega_plot] = meshgrid(linspace(-0.5, 0.5, num_pts_within_roundtrip), ...
    linspace(-1, 1, num_roundtrips - 1) * delta_omega_range);

T = @(x) transmittance(x, M, kappa_TR, varphi_2, varphi_3, u);

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

function y = transmittance(Omegat, M, kappa_TR, varphi_2, varphi_3, u)
U1 = [1, -1; 1, 1]/sqrt(2); 
U2 = [-1i, -1i; 1, -1]/sqrt(2);

y1 = [exp(+1i * kappa_TR * sin(Omegat)), 0;
      0, exp(-1i * kappa_TR * sin(Omegat))];
y2 = [exp(+1i * kappa_TR * sin(M*Omegat+varphi_2)), 0;
      0, exp(-1i * kappa_TR * sin(M*Omegat+varphi_2))];
y3 = [exp(+1i * kappa_TR * (2 + u - cos(Omegat) - cos(M*Omegat+varphi_2) - cos(M^2*Omegat+varphi_3))), 0;
      0, exp(-1i * kappa_TR * (2 + u - cos(Omegat) - cos(M*Omegat+varphi_2) - cos(M^2*Omegat+varphi_3)))];

y = (U1 * y1 * U1' + U2 * y2 * U2' + y3)/3;
end