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

% close all;
clearvars -except dynamics_and_results dynamics_and_results_colourmap;

global legends;
global titles;
global export;
global figures_path;

figures_path = '../figures_output/Figure2/';

legends = true;
titles = false;
export = false;


is_first_correct= true;
normalised = true;


plot_x_pattern(dynamics_and_results,normalised);

plot_accuracy_vs_uncertainty(dynamics_and_results,is_first_correct,normalised)


all = false;
correct =true;
plot_overall_input_to_mc_heatmap(dynamics_and_results, all, correct)

all = false;
correct =false;
plot_overall_input_to_mc_heatmap(dynamics_and_results, all, correct)

function plot_overall_input_to_mc_heatmap(dynamics_and_results_colourmap,all,correct)

global export;
[max_input_0, y_mc_hu_0] = ...
    calculate_mean_activity(dynamics_and_results_colourmap, 0,all,correct);

[max_input_32,y_mc_hu_32] = ...
    calculate_mean_activity(dynamics_and_results_colourmap, 3.2,all,correct);

[max_input_64,y_mc_hu_64] = ...
    calculate_mean_activity(dynamics_and_results_colourmap, 6.4,all,correct);


[max_input_128, y_mc_hu_128] = ...
    calculate_mean_activity(dynamics_and_results_colourmap, 12.8,all,correct);

[max_input_25, y_mc_hu_25] = ...
    calculate_mean_activity(dynamics_and_results_colourmap, 25.6,all,correct);

[max_input_52,y_mc_hu_52] = ...
    calculate_mean_activity(dynamics_and_results_colourmap, 51.2,all,correct);




y=  [max_input_0,max_input_32,max_input_64,max_input_128,max_input_25, max_input_52];

y = y';
z= [y_mc_hu_0,y_mc_hu_32,y_mc_hu_64,y_mc_hu_128,y_mc_hu_25,y_mc_hu_52];
z =z';
n=8000;
x(1:8000)=1;
x(8001:16000)=2;
x(16001:24000)=3;
x(24001:32000)=4;
x(32001:40000)=5;
x(40001:48000)=6;
for i=1:48000
    if(y(i)<0)
        y(i)=0;
    end
end

fullmat = [x;y'];
fullmat = fullmat';
[b2,~] = find(isnan(fullmat));

z(b2,:) = [];
x = x';
x(b2,:) = [];
y(b2,:) = [];

fullmat = [x';y'];
fullmat = fullmat';
[b,indecs]=unique(fullmat,'rows');
z_unique = z(indecs);
x_unique = b(:,1);
y_unique = b(:,2);


n=256;
z_unique= normalise_custom_colour(z_unique,max(z_unique),min(z_unique));
[X,Y] = meshgrid(linspace(min(x_unique),max(x_unique),n), linspace(min(y_unique),max(y_unique),n));

figure8=figure;
hold on;
contourfm(Y,X,griddata(x_unique,y_unique,z_unique,X,Y));

if(~correct)
gcf;    
xlim([1 5]);
xticks([1 2 3 4 5])
xticklabels({'0','3.2','6.4','12.8','25.6'});
end

if(correct)
gcf;
xlim([1 6]);
xticks([1 2 3 4 5 6])
xticklabels({'0','3.2','6.4','12.8','25.6','51.2'});
end

ylim([0 25]);
xlabel(['Evidence quality ',char(949), ' (%)']);
ylabel('Total input to uncertainty population (Hz)');
pubgraph(figure8,28,1,'w');
cb = contourcbar;
cb.XLabel.String  = 'Uncertainty activity level (a.u.)';


label = 'Error responses';
if(correct)
    label='Correct responses';
end

title(label);

cb.XLabel.Position = [1.4400    0.5027         0];
cb.Ticks = [0, 0.99];

set(gcf, 'Position', [240 302 650 545]);

if(export)
export_path = [figures_path  'total_input_uncertainty_all_trials'];
export_fig (export_path, '-nofontswap', '-linecaps','-png', '-transparent','-r300','-q101', '-cmyk','-painters');
savefig(export_path);
end

end

function plot_accuracy_vs_uncertainty(dynamics_and_results,is_first_correct,normalised)
global legends;
global titles;
global export;
global figures_path;
coherences = get_coherence_levels(dynamics_and_results);
n_of_coherences= size(coherences,2);

p_correct = zeros(1,n_of_coherences);
x_correct_max = zeros(1,n_of_coherences);
x_incorrect_max = zeros(1,n_of_coherences);

x_correct_max_std = zeros(1,n_of_coherences);
x_incorrect_max_std = zeros(1,n_of_coherences);
coupled = isfield(dynamics_and_results,'is_first_correct');

is_motor_com = false;
for i=1:n_of_coherences
   coherence= coherences(i);
    
    [noncom_correct_counter, noncom_incorrect_counter,...
        com_correct_counter, com_incorrect_counter, ~, ~, ~] = ...
        calculate_accuracies(dynamics_and_results, coherence, is_first_correct, coupled);
    
    
    
    all_trials = noncom_correct_counter+ noncom_incorrect_counter ...
        +com_correct_counter + com_incorrect_counter;
    [noncom_correct_counter, ~,...
        com_correct_counter, ~, ~, ~, ~] = ...
        calculate_accuracies(dynamics_and_results, coherence, is_first_correct, coupled);
    
    p_correct(i) = (noncom_correct_counter+com_correct_counter)/all_trials;
    [x_correct_max(i), ~, x_correct_max_std(i), ~,...
        ~, ~] = ...
        calculate_x_pattern(dynamics_and_results, coherence, is_motor_com, true);
    
    [x_incorrect_max(i), ~, x_incorrect_max_std(i), ~,...
        ~, ~] = ...
        calculate_x_pattern(dynamics_and_results, coherence, is_motor_com, false);
    
end

confidence_max(1,:)= nanmean([x_correct_max; x_incorrect_max]);
if(normalised)
    confidence_max = normalise_custom(confidence_max, max(confidence_max), min(confidence_max));
end

%%%%%%%%Plot7: Accuracy v.s Uncertainty using peak.

title_string = 'Accuracy v.s uncertainty';
x_string = 'Uncertainty \upsilon';


% [fitObject, gof] = fit(confidence_max',p_correct','poly1');

% 
% Set up fittype and options.
ft = fittype( 'smoothingspline' );
opts = fitoptions( 'Method', 'SmoothingSpline' );
opts.Normalize = 'on';
opts.SmoothingParam = 0.999993297969134;
% 
% % Set up fittype and options.
% ft = fittype( 'smoothingspline' );
% opts = fitoptions( 'Method', 'SmoothingSpline' );
% opts.SmoothingParam = 0.999999666323406;
% 

% Fit model to data.
[fitObject, gof] = fit( confidence_max', p_correct', ft, opts );

%%% extract data from fitted curve
figuretemp = figure;
plot(fitObject);
h = findobj(figuretemp,'Type','line');
x=get(h,'Xdata');
y=get(h,'Ydata');
close;
%%%%%%%


figure7=figure;

plot(x,y,'LineWidth', 4,'Color',[0 0 0]);
hold on;

scatter(confidence_max,p_correct,'o','LineWidth',4,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor', [0 0 0]);
ylim([0.6 1.01])
xlim([0,1.01])
ylabel('Accuracy p');
xlabel(x_string);
pubgraph(figure7,28,4,'w')

if(titles)
    title(title_string);
end

if(export)
export_path = [figures_path  'accuracy_vs_uncertainty'];
export_fig (export_path, '-nofontswap', '-linecaps','-png', '-transparent','-r300','-q101', '-cmyk','-painters');
savefig(export_path);
end
%%%%%%%%

end
function plot_x_pattern(dynamics_and_results,normalised)

global legends;
global titles;
global export;

global figures_path;
coherences = get_coherence_levels(dynamics_and_results);
n_of_coherences= size(coherences,2);

x_correct_max = zeros(1,n_of_coherences);
x_incorrect_max = zeros(1,n_of_coherences);
x_correct_max_std = zeros(1,n_of_coherences);
x_incorrect_max_std = zeros(1,n_of_coherences);

x_correct_integral = zeros(1,n_of_coherences);
x_correct_integral_std = zeros(1,n_of_coherences);
x_incorrect_integral = zeros(1,n_of_coherences);
x_incorrect_integral_std = zeros(1,n_of_coherences);

x_correct_max_lu = zeros(1,n_of_coherences);
x_incorrect_max_lu = zeros(1,n_of_coherences);
x_correct_max_lu_std = zeros(1,n_of_coherences);
x_incorrect_max_lu_std = zeros(1,n_of_coherences);

x_correct_integral_lu = zeros(1,n_of_coherences);
x_correct_integral_lu_std = zeros(1,n_of_coherences);
x_incorrect_integral_lu = zeros(1,n_of_coherences);
x_incorrect_integral_lu_std = zeros(1,n_of_coherences);

is_motor_com= false;

for i=1:n_of_coherences
    coherence = coherences(i);
    [x_correct_max(i), x_correct_max_lu(i), x_correct_max_std(i), x_correct_max_lu_std(i),...
        x_correct_integral(i), x_correct_integral_std(i), x_correct_integral_lu(i), x_correct_integral_lu_std(i)] = ...
        calculate_x_pattern(dynamics_and_results, coherence, is_motor_com, true);
        
    [x_incorrect_max(i), x_incorrect_max_lu(i), x_incorrect_max_std(i), x_incorrect_max_lu_std(i),...
        x_incorrect_integral(i),x_incorrect_integral_std(i), x_incorrect_integral_lu(i), x_incorrect_integral_lu_std(i)] = ...
        calculate_x_pattern(dynamics_and_results, coherence, is_motor_com, false);
end

%%%% HU stuff:
x_max_all = [x_correct_max,x_incorrect_max];
max_bound_peak = max(x_max_all);
min_bound_peak = min(x_max_all);

integral_all = [x_correct_integral,x_incorrect_integral];
max_bound_integral = max(integral_all);
min_bound_integral = min(integral_all);

x_correct_max = normalise_custom(x_correct_max,max_bound_peak,min_bound_peak);
x_correct_max_std = normalise_std(x_correct_max_std,max_bound_peak,min_bound_peak);

x_incorrect_max = normalise_custom(x_incorrect_max,max_bound_peak,min_bound_peak);
x_incorrect_max_std = normalise_std(x_incorrect_max_std,max_bound_peak,min_bound_peak);

x_correct_integral = normalise_custom(x_correct_integral,max_bound_integral,min_bound_integral);
x_correct_integral_std = normalise_std(x_correct_integral_std,max_bound_integral,min_bound_integral);

x_incorrect_integral = normalise_custom(x_incorrect_integral,max_bound_integral,min_bound_integral);
x_incorrect_integral_std = normalise_std(x_incorrect_integral_std,max_bound_integral,min_bound_integral);
%%%%%

%%%LU stuff


x_max_lu_all = [x_correct_max_lu,x_incorrect_max_lu];
max_bound_peak_lu = max(x_max_lu_all);
min_bound_peak_lu = min(x_max_lu_all);

integral_all_lu = [x_correct_integral_lu,x_incorrect_integral_lu];
max_bound_integral_lu = max(integral_all_lu);
min_bound_integral_lu = min(integral_all_lu);

x_correct_max_lu = normalise_custom(x_correct_max_lu,max_bound_peak_lu,min_bound_peak_lu);
x_correct_max_lu_std = normalise_std(x_correct_max_lu_std,max_bound_peak_lu,min_bound_peak_lu);

x_incorrect_max_lu = normalise_custom(x_incorrect_max_lu,max_bound_peak_lu,min_bound_peak_lu);
x_incorrect_max_lu_std = normalise_std(x_incorrect_max_lu_std,max_bound_peak_lu,min_bound_peak_lu);

x_correct_integral_lu = normalise_custom(x_correct_integral_lu,max_bound_integral_lu,min_bound_integral_lu);
x_correct_integral_lu_std = normalise_std(x_correct_integral_lu_std,max_bound_integral_lu,min_bound_integral_lu);

x_incorrect_integral_lu = normalise_custom(x_incorrect_integral_lu,max_bound_integral_lu,min_bound_integral_lu);
x_incorrect_integral_lu_std = normalise_std(x_incorrect_integral_lu_std,max_bound_integral_lu,min_bound_integral_lu);

%%%

%%%%%%%%Plot1: Confidence (X-Pattern) using Peak.


figure1=figure;
axes1 = axes('Parent',figure1);
hold on;
errorbar(x_correct_integral*100, x_correct_integral_std*100, 'LineStyle','-.','LineWidth', 4,'DisplayName',...
    'area, correct','Color',[0         0    0.5430]);

errorbar(x_incorrect_integral*100, x_incorrect_integral_std*100,'LineStyle','-.','LineWidth', 4,'DisplayName',...
    'area, error','Color',[0.9375    0.5000    0.5000]);

errorbar(x_correct_max*100,x_correct_max_std*100, 'LineWidth', 4,'DisplayName',...
    'peak, correct','Color',[0.1172    0.5625    1.0000]);
errorbar(x_incorrect_max*100,x_incorrect_max_std*100, 'LineWidth', 4,'DisplayName',...
    'peak, error','Color',[0.8594    0.0781    0.2344]);

xlim(axes1,[1 6]);
set(axes1,'FontSize',20,'XTickLabel',...
    {'0','3.2','6.4','12.8','25.6','51.2'});

xlabel(['Evidence quality ',char(949), ' (%)']);

title_string = 'Uncertainty v.s Accuracy using Peak';
y_string = 'Peak firing rate of Mean HU Activity (Hz)';

if(normalised)
    title_string = 'Uncertainty v.s Accuracy using Peak Value (Normalised)';
    y_string= 'Uncertainty \upsilon (%)';
end

ylabel(y_string);
pubgraph(figure1,28,4,'w')
ylim([-1,110]);
if(legends)
    legend('show');
end
if(titles)
    title(title_string);
end

if(export)
export_path = [figures_path  'x_pattern'];
export_fig (export_path, '-nofontswap', '-linecaps','-png', '-transparent','-r300','-q101', '-cmyk','-painters');
savefig(export_path);
end
%%%%%%%%

end


function coherences = get_coherence_levels(dynamics_and_results)


coherences = zeros(1,size(dynamics_and_results,1));
for i=1:size(dynamics_and_results)
    coherences(i) = dynamics_and_results(i).coherence_level;
end

coherences = unique(coherences);

return
end


function[max_hu, max_lu, std_hu, std_lu, integral_hu,integral_hu_std, integral_lu,integral_lu_std] = ...
    calculate_x_pattern(dynamics_and_results, coherence, is_motor_com, is_motor_correct)

coherences = get_coherence_levels(dynamics_and_results);

if(~ismember(coherence,coherences))
    error('coherence level not found in data');
end

counter = 0;

trial_length =  dynamics_and_results(1).trial_length;


for i = 1:size(dynamics_and_results,1)
    if(dynamics_and_results(i).coherence_level == coherence &&...
            dynamics_and_results(i).is_motor_com == is_motor_com && ...
            dynamics_and_results(i).is_motor_correct == is_motor_correct ...
            && dynamics_and_results(i).motor_decision_made)
        counter = counter+1;
    end
end

y_mc_hu_gather = zeros(counter,trial_length);
y_mc_lu_gather = zeros(counter,trial_length);

j=1;
for i = 1:size(dynamics_and_results,1)
    if(dynamics_and_results(i).coherence_level == coherence &&...
            dynamics_and_results(i).is_motor_com == is_motor_com && ...
            dynamics_and_results(i).is_motor_correct == is_motor_correct...
            && dynamics_and_results(i).motor_decision_made)
        y_mc_hu_gather(j,:)=dynamics_and_results(i).y_mc_hu;
        y_mc_lu_gather(j,:)=dynamics_and_results(i).y_mc_lu;
        j=j+1;
    end
end


% Max of each trial

maxes_hu_of_each_trial =zeros(1,size(y_mc_hu_gather,1));
maxes_lu_of_each_trial =zeros(1,size(y_mc_hu_gather,1));
for i = 1:size(y_mc_hu_gather,1)
    maxes_hu_of_each_trial(i) = max(y_mc_hu_gather(i,:));
    maxes_lu_of_each_trial(i) = max(y_mc_lu_gather(i,:));
end


%%% calculate x_pattern by peak
max_hu = nanmean(maxes_hu_of_each_trial);
std_hu = nanstd(maxes_hu_of_each_trial)/sqrt(size(maxes_hu_of_each_trial,2));
max_lu = nanmean(maxes_lu_of_each_trial);
std_lu = nanstd(maxes_lu_of_each_trial)/sqrt(size(maxes_lu_of_each_trial,2));



y_mc_hu_mean(1,:) = mean(y_mc_hu_gather);
y_mc_lu_mean(1,:) = mean(y_mc_lu_gather);

%%% calculate x_pattern by integral (trapz)

% Integral of each trial

% integral_hu = trapz(y_mc_hu_mean);
% integral_lu = trapz(y_mc_lu_mean);
integral_hu_of_each_trial =zeros(1,size(y_mc_hu_gather,1));
integral_lu_of_each_trial =zeros(1,size(y_mc_hu_gather,1));
for i = 1:size(y_mc_hu_gather,1)
    integral_hu_of_each_trial(i) = trapz(y_mc_hu_gather(i,:));
    integral_lu_of_each_trial(i) = trapz(y_mc_lu_gather(i,:));
end

integral_hu = nanmean(integral_hu_of_each_trial);
integral_hu_std = nanstd(integral_hu_of_each_trial)/sqrt(size(integral_hu_of_each_trial,2));

integral_lu = nanmean(integral_lu_of_each_trial);
integral_lu_std = nanstd(integral_lu_of_each_trial)/sqrt(size(integral_lu_of_each_trial,2));

return;

end


function [normalised_vector] = normalise_custom(vector,max,min)

normalised_vector = zeros(1,size(vector,2));

for i=1:size(vector,2)
    normalised_vector(i) = (vector(i)-min)/(max-min);
end
return;
end

function [normalised_vector] = normalise_custom_colour(vector,max,min)
vector = vector';
normalised_vector = zeros(1,size(vector,2));

for i=1:size(vector,2)
    normalised_vector(i) = (vector(i)-min)/(max-min);
end
normalised_vector = normalised_vector';
return;
end

function [normalised_std] = normalise_std(std,max,min)

normalised_std= std/(max-min);


return;

end

function[noncom_correct_counter, noncom_incorrect_counter,...
    com_correct_counter, com_incorrect_counter, com_late_correct_counter, com_late_incorrect_counter,nondecision_counter] = ...
    calculate_accuracies(dynamics_and_results, coherence, is_first_correct, coupled)

coherences = get_coherence_levels(dynamics_and_results);

if(~ismember(coherence,coherences))
    error('coherence level not found in data');
end

noncom_correct_counter = 0;
noncom_incorrect_counter = 0;
com_correct_counter = 0;
com_incorrect_counter = 0;

nondecision_counter = 0;

com_late_correct_counter = 0;
com_late_incorrect_counter = 0;

if(~coupled)
for i = 1:size(dynamics_and_results,1)
    if(coherence==dynamics_and_results(i).coherence_level)
        if(dynamics_and_results(i).motor_decision_made)
            if(dynamics_and_results(i).is_motor_correct)
                if(dynamics_and_results(i).is_motor_com)
                    com_correct_counter = com_correct_counter+1;
                    if(dynamics_and_results(i).is_late_com)
                        com_late_correct_counter= com_late_correct_counter +1;
                    end
                else
                    noncom_correct_counter = noncom_correct_counter+1;
                end
            elseif(dynamics_and_results(i).is_motor_com)
                com_incorrect_counter = com_incorrect_counter+1;
                if(dynamics_and_results(i).is_late_com)
                    com_late_incorrect_counter= com_late_incorrect_counter +1;
                end
            else
                noncom_incorrect_counter = noncom_incorrect_counter+1;
            end
        else
            nondecision_counter = nondecision_counter+1;
        end
    end
end
return;
end

for i = 1:size(dynamics_and_results,1)
    if(coherence==dynamics_and_results(i).coherence_level)
        if(is_first_correct==dynamics_and_results(i).is_first_correct)
            if(dynamics_and_results(i).motor_decision_made)
                if(dynamics_and_results(i).is_motor_correct)
                    if(dynamics_and_results(i).is_motor_com)
                        com_correct_counter = com_correct_counter+1;
                        if(dynamics_and_results(i).is_late_com)
                            com_late_correct_counter= com_late_correct_counter +1;
                        end
                    else
                        noncom_correct_counter = noncom_correct_counter+1;
                    end
                elseif(dynamics_and_results(i).is_motor_com)
                    com_incorrect_counter = com_incorrect_counter+1;
                    if(dynamics_and_results(i).is_late_com)
                        com_late_incorrect_counter= com_late_incorrect_counter +1;
                    end
                else
                    noncom_incorrect_counter = noncom_incorrect_counter+1;
                end
            else
                nondecision_counter = nondecision_counter+1;
            end
        end
    end
end
return;


end

function[mc_input_max, mc_input_max_std,uncertainty, uncertainty_std] = ...
    calculate_max_input(dynamics_and_results, coherence)

coherences = get_coherence_levels(dynamics_and_results);

if(~ismember(coherence,coherences))
    error('coherence level not found in data');
end

counter = 0;

trial_length =  dynamics_and_results(1).trial_length;


for i = 1:size(dynamics_and_results,1)
    if(dynamics_and_results(i).coherence_level == coherence)
        counter = counter+1;
    end
end

max_input_gather = zeros(counter,trial_length);
uncertainty_gather = zeros(counter,trial_length);

j=1;
for i = 1:size(dynamics_and_results,1)
    if(dynamics_and_results(i).coherence_level == coherence)
        max_input_gather(j,:)=dynamics_and_results(i).overall_mc_input;
        uncertainty_gather(j,:)=dynamics_and_results(i).y_mc_hu;
        j=j+1;
    end
end


% Max of each trial

maxes_of_mc_input =zeros(1,size(max_input_gather,1));
maxes_of_uncertainty =zeros(1,size(uncertainty_gather,1));
for i = 1:size(max_input_gather,1)
    maxes_of_mc_input(i) = max(max_input_gather(i,:));
    maxes_of_uncertainty(i) = max(uncertainty_gather(i,:));
    
end


%%% calculate x_pattern by peak
mc_input_max = nanmean(maxes_of_mc_input);
mc_input_max_std = nanstd(maxes_of_mc_input)/sqrt(size(maxes_of_mc_input,2));

%%% calculate x_pattern by peak
uncertainty = nanmean(maxes_of_uncertainty);
uncertainty_std = nanstd(maxes_of_uncertainty)/sqrt(size(maxes_of_uncertainty,2));

return;

end


function[overall_mc_input_mean,y_mc_hu_mean] = ...
    calculate_mean_activity(dynamics_and_results, coherence,all,correct)

coherences = get_coherence_levels(dynamics_and_results);

if(~ismember(coherence,coherences))
    error('coherence level not found in data');
end

counter = 0;
trial_length =  dynamics_and_results(1).trial_length;

%get number of each type, before initialising vectors.
if(all)
    
    for i = 1:size(dynamics_and_results,1)
        if(dynamics_and_results(i).coherence_level == ...
                coherence && dynamics_and_results(i).motor_decision_made)
            counter = counter+1;
        end
    end
    
    max_input_gather = zeros(counter,trial_length);
    y_mc_hu_gather = zeros(counter,trial_length);
    
    j=1;
    for i = 1:size(dynamics_and_results,1)
        if(dynamics_and_results(i).coherence_level == ...
                coherence && dynamics_and_results(i).motor_decision_made)
            max_input_gather(j,:)=dynamics_and_results(i).overall_mc_input;
            y_mc_hu_gather(j,:)=dynamics_and_results(i).y_mc_hu;
            j=j+1;
        end
    end
    
else
    for i = 1:size(dynamics_and_results,1)
        if(dynamics_and_results(i).coherence_level == coherence && ...
                dynamics_and_results(i).is_motor_correct == correct &&...
                dynamics_and_results(i).motor_decision_made)
            counter = counter+1;
        end
    end
    
    max_input_gather = zeros(counter,trial_length);
    y_mc_hu_gather = zeros(counter,trial_length);
    overall_mc_input_mean = zeros(1, trial_length);
    y_mc_hu_mean = zeros(1, trial_length);
    j=1;
    for i = 1:size(dynamics_and_results,1)
        if(dynamics_and_results(i).coherence_level == coherence && ...
                dynamics_and_results(i).is_motor_correct == correct &&...
                dynamics_and_results(i).motor_decision_made)
            max_input_gather(j,:)=dynamics_and_results(i).overall_mc_input;
            y_mc_hu_gather(j,:)=dynamics_and_results(i).y_mc_hu;
            j=j+1;
        end
    end
    
end

% Non CoM Correct Mean
% overall_mc_input_mean(1,:) = mean(max_input_gather);
% y_mc_hu_mean(1,:) = mean(y_mc_hu_gather);

if(~isempty(max_input_gather))
    overall_mc_input_mean = max_input_gather(1,:);
    y_mc_hu_mean = y_mc_hu_gather(1,:);
end


return;
end


