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

global legends;
global titles;
global export;
global figures_path;

figures_path = '../figures_output/Figure4/';

legends = true;
titles = false;
export = false;

normalised=true;
is_first_correct= true;

plot_accuracies_and_pcom(dynamics_and_results,is_first_correct,normalised);

function plot_accuracies_and_pcom(dynamics_and_results, is_first_correct,normalised)
global legends;
global titles;
global export;
global figures_path;
coherences = get_coherence_levels(dynamics_and_results);
n_of_coherences= size(coherences,2);

p_com = zeros(1,n_of_coherences);
p_com_correct = zeros(1,n_of_coherences);
p_com_incorrect = zeros(1,n_of_coherences);

com_correct = zeros(1,n_of_coherences);
com_incorrect = zeros(1,n_of_coherences);

p_late_com_correct = zeros(1,n_of_coherences);
p_late_com_incorrect = zeros(1,n_of_coherences);

late_com_correct = zeros(1,n_of_coherences);
late_com_incorrect = zeros(1,n_of_coherences);

p_correct = zeros(1,n_of_coherences);
p_correct_noncom = zeros(1,n_of_coherences);
p_nondecision = zeros(1,n_of_coherences);
initiation_time_noncom_correct_mean_gather = zeros(1,n_of_coherences);
initiation_time_noncom_incorrect_mean_gather = zeros(1,n_of_coherences);
initiation_time_com_correct_mean_gather = zeros(1,n_of_coherences);
initiation_time_com_incorrect_mean_gather = zeros(1,n_of_coherences);

initiation_time_noncom_correct_std_gather = zeros(1,n_of_coherences);
initiation_time_noncom_incorrect_std_gather = zeros(1,n_of_coherences);
initiation_time_com_correct_std_gather = zeros(1,n_of_coherences);
initiation_time_com_incorrect_std_gather = zeros(1,n_of_coherences);

coupled = isfield(dynamics_and_results,'is_first_correct');

for i=1:n_of_coherences
    coherence= coherences(i);
    
    [noncom_correct_counter, noncom_incorrect_counter,...
        com_correct_counter, com_incorrect_counter, com_late_correct_counter, com_late_incorrect_counter, nondecision_counter] = ...
        calculate_accuracies(dynamics_and_results, coherence, is_first_correct, coupled);
    
    
    
    all_trials = noncom_correct_counter+ noncom_incorrect_counter ...
        +com_correct_counter + com_incorrect_counter;
    p_com_incorrect(i) = com_incorrect_counter/all_trials;
    p_com_correct(i) = com_correct_counter/all_trials;
    com_correct(i) = com_correct_counter;
    com_incorrect(i) = com_incorrect_counter;
    p_late_com_correct(i) = com_late_correct_counter/all_trials;
    p_late_com_incorrect(i) = com_late_incorrect_counter/all_trials;
    late_com_correct(i) = com_late_correct_counter;
    late_com_incorrect(i) = com_late_incorrect_counter;
    p_com(i)=(com_correct_counter+com_incorrect_counter)/all_trials;
    p_correct(i) = (noncom_correct_counter+com_correct_counter)/all_trials;
    p_correct_noncom(i) = (noncom_correct_counter)/all_trials;
    p_nondecision(i)= nondecision_counter/all_trials;
    [initiation_time_noncom_correct_mean_gather(i),initiation_time_noncom_correct_std_gather(i),...
        ~,~] = ...
        calculate_mean_times(dynamics_and_results, coherence, ...
        true, false, is_first_correct, coupled);
    
    [initiation_time_noncom_incorrect_mean_gather(i),initiation_time_noncom_incorrect_std_gather(i),...
        ~,~] = ...
        calculate_mean_times(dynamics_and_results, coherence, ...
        false, false, is_first_correct, coupled);
    
    [initiation_time_com_correct_mean_gather(i),initiation_time_com_correct_std_gather(i),...
        ~,~] = ...
        calculate_mean_times(dynamics_and_results, coherence, ...
        true, true, is_first_correct, coupled);
    
    [initiation_time_com_incorrect_mean_gather(i),initiation_time_com_incorrect_std_gather(i),...
        ~,~] = ...
        calculate_mean_times(dynamics_and_results, coherence, ...
        false, true, is_first_correct, coupled);
    
end

early_com_correct = com_correct - late_com_correct;
early_com_incorrect = com_incorrect - late_com_incorrect;

p_early_com_correct = early_com_correct./all_trials;
p_early_com_incorrect = early_com_incorrect./all_trials;

% p_com = p_early_com_correct + p_early_com_incorrect;
% p_com_correct = p_early_com_correct;
% p_com_incorrect = p_early_com_incorrect;
%%%%%%%%Plot 1: PCoM related results
figure1=figure;
axes1 = axes('Parent',figure1);
plot(p_com, 'LineWidth', 4,'DisplayName',...
    'all','Color',[0 0 0]);
hold on;
plot(p_com_correct, 'LineWidth', 4,'DisplayName',...
    'correct','Color',[0 1 0]);
plot(p_com_incorrect, 'LineWidth', 4,'DisplayName',...
    'error','Color',[1 0 0]);
xlim(axes1,[1 6]);

set(axes1,'FontSize',20,'XTickLabel',...
    {'0','3.2','6.4','12.8','25.6','51.2'});
xlabel(['Evidence quality ',char(949), ' (%)']);
ylabel('Probability of change-of-mind');
pubgraph(figure1,28,4,'w')

title_string= 'Probability of change-of-mind';
if(legends)
    legend('show');
end
if(titles)
    title(title_string);
end

if(export)
export_path = [figures_path  'p_com_errorvscorrect'];
export_fig (export_path, '-nofontswap', '-linecaps','-png', '-transparent','-r300','-q101', '-cmyk','-painters');
savefig(export_path);
end
set(gcf, 'Position', [779 304 650 500])
%%%%%%%%

%%%%%%%%Plot 5: Late PCoM related results
figure5=figure;
axes1 = axes('Parent',figure5);
hold on;
plot(p_late_com_correct*100, 'LineWidth', 4,'DisplayName',...
    'late, correct','Color',[0.1250    0.6953    0.6641]);

plot(p_late_com_incorrect*100, 'LineWidth', 4,'DisplayName',...
    'late, error','Color',[0.9102    0.5859    0.4766]);

plot(p_early_com_correct*100, 'LineStyle','-.','LineWidth', 4,'DisplayName',...
    'early, correct','Color',[0.4844    0.9875         0]);

plot(p_early_com_incorrect*100,'LineStyle','-.', 'LineWidth', 4,'DisplayName',...
    'early, error','Color',[0.5430         0         0]);


xlim(axes1,[1 6]);

set(axes1,'FontSize',20,'XTickLabel',...
    {'0','3.2','6.4','12.8','25.6','51.2'});
xlabel(['Evidence quality ',char(949), ' (%)']);
ylabel('Probability of change-of-mind');
pubgraph(figure5,28,4,'w')
title_string = 'probability of late change';
if(legends)
    legend('show');
end
if(titles)
    title(title_string);
end

if(export)
export_path = [figures_path  'p_com_late_vs_early'];
export_fig (export_path, '-nofontswap', '-linecaps','-png', '-transparent','-r300','-q101', '-cmyk','-painters');
savefig(export_path);
end
set(gcf, 'Position', [779 304 650 500])

%%%%%%%% Plot 8: Response time

all_times = [initiation_time_noncom_correct_mean_gather,initiation_time_noncom_incorrect_mean_gather,initiation_time_com_correct_mean_gather,initiation_time_com_incorrect_mean_gather];


normalised_initiation_time_noncom_correct = normalise_custom(initiation_time_noncom_correct_mean_gather, max(all_times),min(all_times));
normalised_initiation_time_noncom_correct_std = normalise_std(initiation_time_noncom_correct_std_gather, max(all_times),min(all_times));

normalised_initiation_time_noncom_incorrect = normalise_custom(initiation_time_noncom_incorrect_mean_gather, max(all_times),min(all_times));
normalised_initiation_time_noncom_incorrect_std = normalise_std(initiation_time_noncom_incorrect_std_gather, max(all_times),min(all_times));


normalised_initiation_time_com_correct_mean_gather = normalise_custom(initiation_time_com_correct_mean_gather, max(all_times),min(all_times));
normalised_initiation_time_com_correct_std_gather = normalise_std(initiation_time_com_correct_std_gather, max(all_times),min(all_times));

normalised_initiation_time_com_incorrect_mean_gather = normalise_custom(initiation_time_com_incorrect_mean_gather, max(all_times),min(all_times));
normalised_initiation_time_com_incorrect_std_gather = normalise_std(initiation_time_com_incorrect_std_gather, max(all_times),min(all_times));

figure8=figure;
axes1 = axes('Parent',figure8);
hold on;
errorbar(normalised_initiation_time_noncom_correct,...
    normalised_initiation_time_noncom_correct_std...
    ,'LineWidth', 4,'DisplayName',...
    'non-com, correct','Color',[0 0 1]);
errorbar(normalised_initiation_time_noncom_incorrect,...
    normalised_initiation_time_noncom_incorrect_std...
    ,'LineWidth', 4,'DisplayName',...
    'non-com, error','Color',[1 0 0]);
errorbar(normalised_initiation_time_com_correct_mean_gather,...
    normalised_initiation_time_com_correct_std_gather...
    ,'LineWidth', 6,'DisplayName',...
    'com, correct','Color',[0  0.3906  0],'LineStyle',':');
errorbar(normalised_initiation_time_com_incorrect_mean_gather,...
    normalised_initiation_time_com_incorrect_std_gather...
    ,'LineWidth', 2,'DisplayName',...
    'com, error','Color',[0.5430  0   0],'LineStyle','-');

xlim(axes1,[1 6]);
ylim([0 1]);
set(axes1,'FontSize',20,'XTickLabel',...
    {'0','3.2','6.4','12.8','25.6','51.2'});

xlabel(['Evidence quality ',char(949), ' (%)']);
ylabel('Normalized response time');
pubgraph(figure8,28,4,'w')
title_string = 'Response time v.s coherence level';
if(legends)
    legend('show');
end
if(titles)
    title(title_string);
end

if(export)
export_path = [figures_path  'reactiontime_noncomvscom'];
export_fig (export_path, '-nofontswap', '-linecaps','-png', '-transparent','-r300','-q101', '-cmyk','-painters');
savefig(export_path);
end

set(gcf, 'Position', [779 304 650 500])

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
function[initiation_time_mean, initiation_time_std, ...
    response_time_mean, response_time_std] = ...
    calculate_mean_times(dynamics_and_results, coherence, ...
    is_motor_correct, is_motor_com, is_first_correct, coupled)

coherences = get_coherence_levels(dynamics_and_results);

if(~ismember(coherence,coherences))
    error('coherence level not found in data');
end

counter = 0;

if(~coupled)
    for i = 1:size(dynamics_and_results,1)
        if(dynamics_and_results(i).coherence_level ==...
                coherence && dynamics_and_results(i).is_motor_correct ==...
                is_motor_correct && dynamics_and_results(i).is_motor_com == ...
                is_motor_com && dynamics_and_results(i).motor_decision_made)
            counter = counter+1;
        end
    end
else
    for i = 1:size(dynamics_and_results,1)
        if(dynamics_and_results(i).coherence_level ==...
                coherence && dynamics_and_results(i).is_motor_correct ==...
                is_motor_correct && dynamics_and_results(i).is_motor_com == ...
                is_motor_com && dynamics_and_results(i).is_first_correct == ...
                is_first_correct && dynamics_and_results(i).motor_decision_made)
            counter = counter+1;
        end
    end
end

if(counter<5)
    initiation_time_mean = NaN;
    initiation_time_std = NaN;
    response_time_mean = NaN;
    response_time_std = NaN;
    return;
end

initiation_time_gather = zeros(1,counter);
response_time_gather = zeros(1,counter);

j=1;
if(~coupled)
    for i = 1:size(dynamics_and_results,1)
        if(dynamics_and_results(i).coherence_level == ...
                coherence && dynamics_and_results(i).is_motor_correct ==...
                is_motor_correct && dynamics_and_results(i).is_motor_com ==...
                is_motor_com && dynamics_and_results(i).motor_decision_made)
            initiation_time_gather(j)=dynamics_and_results(i).response_time -dynamics_and_results(i).stim_onset ;
            response_time_gather(j)=dynamics_and_results(i).response_time + dynamics_and_results(i).movement_duration;
            j=j+1;
        end
    end
else
    for i = 1:size(dynamics_and_results,1)
        if(dynamics_and_results(i).coherence_level == ...
                coherence && dynamics_and_results(i).is_motor_correct ==...
                is_motor_correct && dynamics_and_results(i).is_motor_com ==...
                is_motor_com && dynamics_and_results(i).is_first_correct ==...
                is_first_correct && dynamics_and_results(i).motor_decision_made)
            initiation_time_gather(j)=dynamics_and_results(i).response_time -dynamics_and_results(i).stim_onset ;
            response_time_gather(j)=dynamics_and_results(i).response_time + dynamics_and_results(i).movement_duration;
            j=j+1;
        end
    end
end

initiation_time_mean = mean(initiation_time_gather);
initiation_time_std = std(initiation_time_gather)/sqrt(size(initiation_time_gather,2));
response_time_mean = mean(response_time_gather);
response_time_std = std(response_time_gather)/sqrt(size(response_time_gather,2));


return;
end


function [normalised_vector] = normalise_custom(vector,max,min)

normalised_vector = zeros(1,size(vector,2));

for i=1:size(vector,2)
    normalised_vector(i) = (vector(i)-min)/(max-min);
end
return;
end

function [normalised_std] = normalise_std(std,max,min)

normalised_std= std/(max-min);


return;

end