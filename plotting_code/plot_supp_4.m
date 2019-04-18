% Copyright 2019 Nadim Atiya
% 
% uncertainty_com_modelling is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% uncertainty_com_modelling is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with uncertainty_com_modelling.  If not, see <https://www.gnu.org/licenses/>.

close all;

plot_p_com_colormap(analysis_results)
plot_p_correct_colormap(analysis_results)
plot_rt_colormap(analysis_results)

function plot_p_com_colormap(analysis_results)
n=512;
j = extractfield(analysis_results,'j_mc_value');
mu = extractfield(analysis_results, 'constant_input_value');
[X,Y] = meshgrid(linspace(min(j),max(j),n), linspace(min(mu),max(mu),n));

z_unique = extractfield(analysis_results,'p_com');

figure8=figure;
hold on;
contourfm(Y,X,griddata(j,mu,z_unique,X,Y));

ylim([min(mu) max(mu)]);
xlim([min(j) max(j)]);
xlabel('J_{mc0} (nA Hz^{-1})');
ylabel('\mu (Hz)');
cb = contourcbar;
cb.XLabel.String  = 'P(CoM)';
cb.Ticks = [0, 0.02,0.07,0.09];
cb.XLabel.Position= [1.9088 0.0486 0]
pubgraph(figure8,28,1,'w');


set(gcf, 'Position', [291 457 535 420])
end

function plot_p_correct_colormap(analysis_results)
n=256;
j = extractfield(analysis_results,'j_mc_value');
mu = extractfield(analysis_results, 'constant_input_value');
[X,Y] = meshgrid(linspace(min(j),max(j),n), linspace(min(mu),max(mu),n));

z_unique = extractfield(analysis_results,'p_correct');

figure8=figure;
hold on;
contourfm(Y,X,griddata(j,mu,z_unique,X,Y));

ylim([min(mu) max(mu)]);
xlim([min(j) max(j)]);
xlabel('J_{mc0} (nA Hz^{-1})');
ylabel('\mu (Hz)');



pubgraph(figure8,28,1,'w');
cb = contourcbar;
cb.XLabel.String  = 'P(Correct)';
cb.Ticks = [0.76, 0.78, 0.86, 0.88];
cb.XLabel.Position = [1.9713 0.8220 0];
set(gcf, 'Position', [291 457 535 420])
end

function plot_rt_colormap(analysis_results)
n=256;
j = extractfield(analysis_results,'j_mc_value');
mu = extractfield(analysis_results, 'constant_input_value');
[X,Y] = meshgrid(linspace(min(j),max(j),n), linspace(min(mu),max(mu),n));

z_unique = extractfield(analysis_results,'mean_rt');

figure8=figure;
hold on;
contourfm(Y,X,griddata(j,mu,z_unique,X,Y));

ylim([min(mu) max(mu)]);
xlim([min(j) max(j)]);
xlabel('J_{mc0} (nA Hz^{-1})');
ylabel('\mu (Hz)');

pubgraph(figure8,28,1,'w');
cb = contourcbar;
cb.Ticks = [600,700,1150,1250];

cb.XLabel.Position = [2.0594  928.6176         0];
cb.XLabel.String  = 'RT (ms)';

set(gcf, 'Position', [291 457 535 420])
end