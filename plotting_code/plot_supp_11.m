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
n=256;
j = extractfield(analysis_results,'inhibitory_onset');
mu = extractfield(analysis_results, 'uncertainty_onset');
[X,Y] = meshgrid(linspace(min(j),max(j),n), linspace(min(mu),max(mu),n));

z_unique = extractfield(analysis_results,'p_com');

figure8=figure;
hold on;
clev = [0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15,0.16,0.17];
contourfcmap(X,Y,griddata(j,mu,z_unique,X,Y),clev,jet(14), ...
    'method', 'recolor');
ylim([min(mu) max(mu)]);
xlim([min(j) max(j)]);
xlabel('t_{INH}(ms)');
ylabel('t_{U}(ms)');
cb = contourcbar;
cb.Ticks = [0.04, 0.06,0.14, 0.16];
cb.XLabel.String  = 'P(CoM)';
cb.XLabel.Position = [1.4088    0.1039         0];

pubgraph(figure8,28,1,'w');
set(gcf, 'Position', [79 422 560 420])

end

function plot_p_correct_colormap(analysis_results)
n=256;
j = extractfield(analysis_results,'inhibitory_onset');
mu = extractfield(analysis_results, 'uncertainty_onset');
[X,Y] = meshgrid(linspace(min(j),max(j),n), linspace(min(mu),max(mu),n));

z_unique = extractfield(analysis_results,'p_correct');
clev = [0.57,0.6,0.63,0.66,0.69,0.72,0.75,0.78,0.81];
figure8=figure;
hold on;
contourfcmap(X,Y,griddata(j,mu,z_unique,X,Y),clev,jet(8), ...
    'method', 'recolor')

ylim([min(mu) max(mu)]);
xlim([min(j) max(j)]);
xlabel('t_{INH}(ms)');
ylabel('t_{U}(ms)');
cb = contourcbar;
cb.XLabel.String  = 'P(Correct)';
cb.Ticks = [0.58,0.61 0.75,0.78];
cb.XLabel.Position = [1.9713    0.6811         0];

pubgraph(figure8,28,1,'w');
set(gcf, 'Position', [79 422 560 420])
end

function plot_rt_colormap(analysis_results)
n=256;
j = extractfield(analysis_results,'inhibitory_onset');
mu = extractfield(analysis_results, 'uncertainty_onset');
[X,Y] = meshgrid(linspace(min(j),max(j),n), linspace(min(mu),max(mu),n));

z_unique = extractfield(analysis_results,'mean_rt');
clev = [50,100,150,200,250,300,350,400,450,500,550,600,650,700];

figure8=figure;
hold on;
contourfcmap(X,Y,griddata(j,mu,z_unique,X,Y),clev,jet(13), ...
    'method', 'recolor')

ylim([min(mu) max(mu)]);
xlim([min(j) max(j)]);
xlabel('t_{INH}(ms)');
ylabel('t_{U}(ms)');
cb = contourcbar;
cb.XLabel.Position = [1.8463  350.5200         0];
cb.Ticks = [100, 200, 500, 600];
cb.XLabel.String  = 'RT (ms)';

pubgraph(figure8,28,1,'w');
set(gcf, 'Position', [79 422 560 420])
end