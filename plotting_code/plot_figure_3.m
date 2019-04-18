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

% 
% %call functions to test:
close all;
clearvars -except dynamics_and_results dynamics_and_results_2nddecision;

global legends;
global titles;
global export;
global figures_path;

figures_path = '../figures_output/Figure3/';

legends = true;
titles = false;
export = false;


is_first_correct= true;
normalised = true;

plot_kiani_response_time(dynamics_and_results_2nddecision)

function plot_kiani_response_time(dynamics_and_results_2nddecision)

global legends;
global titles;
global export;
global figures_path;
coherences = get_coherence_levels(dynamics_and_results_2nddecision);
n_of_coherences= size(coherences,2);

initiation_time_noncom_correct_mean_gather_2nd = zeros(1,n_of_coherences);
initiation_time_noncom_incorrect_mean_gather_2nd = zeros(1,n_of_coherences);
initiation_time_com_correct_mean_gather_2nd = zeros(1,n_of_coherences);
initiation_time_com_incorrect_mean_gather_2nd = zeros(1,n_of_coherences);

initiation_time_noncom_correct_std_gather_2nd = zeros(1,n_of_coherences);
initiation_time_noncom_incorrect_std_gather_2nd = zeros(1,n_of_coherences);
initiation_time_com_correct_std_gather_2nd = zeros(1,n_of_coherences);
initiation_time_com_incorrect_std_gather_2nd = zeros(1,n_of_coherences);


initiation_time_noncom_correct_mean_gather_2nd_incorrect = zeros(1,n_of_coherences);
initiation_time_noncom_incorrect_mean_gather_2nd_incorrect = zeros(1,n_of_coherences);
initiation_time_com_correct_mean_gather_2nd_incorrect = zeros(1,n_of_coherences);
initiation_time_com_incorrect_mean_gather_2nd_incorrect = zeros(1,n_of_coherences);

initiation_time_noncom_correct_std_gather_2nd_incorrect = zeros(1,n_of_coherences);
initiation_time_noncom_incorrect_std_gather_2nd_incorrect = zeros(1,n_of_coherences);
initiation_time_com_correct_std_gather_2nd_incorrect = zeros(1,n_of_coherences);
initiation_time_com_incorrect_std_gather_2nd_incorrect = zeros(1,n_of_coherences);


response_time_noncom_correct_mean_gather_2nd = zeros(1,n_of_coherences);
response_time_noncom_incorrect_mean_gather_2nd = zeros(1,n_of_coherences);
response_time_com_correct_mean_gather_2nd = zeros(1,n_of_coherences);
response_time_com_incorrect_mean_gather_2nd = zeros(1,n_of_coherences);

response_time_noncom_correct_std_gather_2nd = zeros(1,n_of_coherences);
response_time_noncom_incorrect_std_gather_2nd = zeros(1,n_of_coherences);
response_time_com_correct_std_gather_2nd = zeros(1,n_of_coherences);
response_time_com_incorrect_std_gather_2nd = zeros(1,n_of_coherences);


response_time_noncom_correct_mean_gather_2nd_incorrect = zeros(1,n_of_coherences);
response_time_noncom_incorrect_mean_gather_2nd_incorrect = zeros(1,n_of_coherences);
response_time_com_correct_mean_gather_2nd_incorrect = zeros(1,n_of_coherences);
response_time_com_incorrect_mean_gather_2nd_incorrect = zeros(1,n_of_coherences);

response_time_noncom_correct_std_gather_2nd_incorrect = zeros(1,n_of_coherences);
response_time_noncom_incorrect_std_gather_2nd_incorrect = zeros(1,n_of_coherences);
response_time_com_correct_std_gather_2nd_incorrect = zeros(1,n_of_coherences);
response_time_com_incorrect_std_gather_2nd_incorrect = zeros(1,n_of_coherences);


for i=1:n_of_coherences  
    coherence= coherences(i);
    coupled = true;
    is_first_correct = true;
    
    [initiation_time_noncom_correct_mean_gather_2nd(i),initiation_time_noncom_correct_std_gather_2nd(i),...
        response_time_noncom_correct_mean_gather_2nd(i),response_time_noncom_correct_std_gather_2nd(i)] = ...
        calculate_mean_times(dynamics_and_results_2nddecision, coherence, ...
        true, false, is_first_correct, coupled);
    
    [initiation_time_noncom_incorrect_mean_gather_2nd(i),initiation_time_noncom_incorrect_std_gather_2nd(i),...
        response_time_noncom_incorrect_mean_gather_2nd(i),response_time_noncom_incorrect_std_gather_2nd(i)] = ...
        calculate_mean_times(dynamics_and_results_2nddecision, coherence, ...
        false, false, is_first_correct, coupled);
    
    [initiation_time_com_correct_mean_gather_2nd(i),initiation_time_com_correct_std_gather_2nd(i),...
        response_time_com_correct_mean_gather_2nd(i),response_time_com_correct_std_gather_2nd(i)] = ...
        calculate_mean_times(dynamics_and_results_2nddecision, coherence, ...
        true, true, is_first_correct, coupled);
    
    [initiation_time_com_incorrect_mean_gather_2nd(i),initiation_time_com_incorrect_std_gather_2nd(i),...
        response_time_com_incorrect_mean_gather_2nd(i),response_time_com_incorrect_std_gather_2nd(i)] = ...
        calculate_mean_times(dynamics_and_results_2nddecision, coherence, ...
        false, true, is_first_correct, coupled);
   
    is_first_correct = false;
    
        [initiation_time_noncom_correct_mean_gather_2nd_incorrect(i),initiation_time_noncom_correct_std_gather_2nd_incorrect(i),...
        response_time_noncom_correct_mean_gather_2nd_incorrect(i),response_time_noncom_correct_std_gather_2nd_incorrect(i)] = ...
        calculate_mean_times(dynamics_and_results_2nddecision, coherence, ...
        true, false, is_first_correct, coupled);
    
    [initiation_time_noncom_incorrect_mean_gather_2nd_incorrect(i),initiation_time_noncom_incorrect_std_gather_2nd_incorrect(i),...
        response_time_noncom_incorrect_mean_gather_2nd_incorrect(i),response_time_noncom_incorrect_std_gather_2nd_incorrect(i)] = ...
        calculate_mean_times(dynamics_and_results_2nddecision, coherence, ...
        false, false, is_first_correct, coupled);
    
    [initiation_time_com_correct_mean_gather_2nd_incorrect(i),initiation_time_com_correct_std_gather_2nd_incorrect(i),...
        response_time_com_correct_mean_gather_2nd_incorrect(i),response_time_com_correct_std_gather_2nd_incorrect(i)] = ...
        calculate_mean_times(dynamics_and_results_2nddecision, coherence, ...
        true, true, is_first_correct, coupled);
    
    [initiation_time_com_incorrect_mean_gather_2nd_incorrect(i),initiation_time_com_incorrect_std_gather_2nd_incorrect(i),...
        response_time_com_incorrect_mean_gather_2nd_incorrect(i),response_time_com_incorrect_std_gather_2nd_incorrect(i)] = ...
        calculate_mean_times(dynamics_and_results_2nddecision, coherence, ...
        false, true, is_first_correct, coupled);
end

total_noncom_2nd(1,:) = nanmean([initiation_time_noncom_correct_mean_gather_2nd;  initiation_time_noncom_incorrect_mean_gather_2nd]);
total_noncom_std_2nd(1,:) = nanmean([initiation_time_noncom_correct_std_gather_2nd;  initiation_time_noncom_incorrect_std_gather_2nd]);

total_noncom_response_2nd(1,:) = nanmean([response_time_noncom_correct_mean_gather_2nd;  response_time_noncom_incorrect_mean_gather_2nd]);
total_noncom_response_std_2nd(1,:) = nanmean([response_time_noncom_correct_std_gather_2nd;  response_time_noncom_incorrect_std_gather_2nd]);


total_noncom_2nd_incorrect(1,:) = nanmean([initiation_time_noncom_correct_mean_gather_2nd_incorrect;  initiation_time_noncom_incorrect_mean_gather_2nd_incorrect]);
total_noncom_std_2nd_incorrect(1,:) = nanmean([initiation_time_noncom_correct_std_gather_2nd_incorrect;  initiation_time_noncom_incorrect_std_gather_2nd_incorrect]);

total_noncom_response_2nd_incorrect(1,:) = nanmean([response_time_noncom_correct_mean_gather_2nd_incorrect;  response_time_noncom_incorrect_mean_gather_2nd_incorrect]);
total_noncom_response_std_2nd_incorrect(1,:) = nanmean([response_time_noncom_correct_std_gather_2nd_incorrect;  response_time_noncom_incorrect_std_gather_2nd_incorrect]);

%%%%%%%%Plot1: Overall  time (1st correct dec vs 1st incorrect)
figure1=figure;
axes1 = axes('Parent',figure1);
hold on;


total = [total_noncom_response_2nd, total_noncom_response_2nd_incorrect];

total_noncom_response_2nd_normalised = normalise_custom(total_noncom_response_2nd, max(total),min(total));
total_noncom_response_std_2nd_normalised = normalise_std(total_noncom_response_std_2nd, max(total),min(total));
total_noncom_response_2nd_incorrect_normalised = normalise_custom(total_noncom_response_2nd_incorrect, max(total),min(total));
total_noncom_response_std_2nd_incorrect_normalised = normalise_std(total_noncom_response_std_2nd_incorrect, max(total),min(total));


errorbar(total_noncom_response_2nd_normalised,...
    total_noncom_response_std_2nd_normalised...
    ,'LineWidth', 4,'DisplayName',...
    '1^{st} decision corect','Color',[0.0977    0.0977    0.4375]);
errorbar(total_noncom_response_2nd_incorrect_normalised,...
    total_noncom_response_std_2nd_incorrect_normalised...
    ,'LineWidth', 4,'DisplayName',...
    '1^{st} decision error','Color',[0.6953    0.1328    0.1328]);


xlim(axes1,[1 6]);
ylim([0 1])
set(axes1,'FontSize',20,'XTickLabel',...     {'0','3.2','6.4','12.8','25.6','51.2'});
    {'0','3.2','6.4','12.8','25.6','51.2'});


xlabel(['Evidence quality ',char(949), ' (%)']);
ylabel('Normalized response time of 2^{nd} decision');
pubgraph(figure1,28,4,'w')
title_string = 'Response times of 2^{nd} decision';
if(legends)
    legend('show');
end
if(titles)
    title(title_string);
end

set(gcf, 'Position', [161 291 640 680])

if(export)
export_path = [figures_path  'kiani_overall'];
export_fig (export_path, '-nofontswap', '-linecaps','-png', '-transparent','-r300','-q101', '-cmyk','-painters');
savefig(export_path);
end

%%%%%%%%

%%%%%%%%Plot2: Initiation time (detailed)
all_responses = [response_time_noncom_correct_mean_gather_2nd,response_time_noncom_incorrect_mean_gather_2nd,response_time_noncom_correct_mean_gather_2nd_incorrect,response_time_noncom_incorrect_mean_gather_2nd_incorrect];
response_time_noncom_correct_mean_gather_2nd_normalised = normalise_custom(response_time_noncom_correct_mean_gather_2nd, max(all_responses),min(all_responses));
response_time_noncom_correct_std_gather_2nd_normalised = normalise_std(response_time_noncom_correct_std_gather_2nd, max(all_responses),min(all_responses));

response_time_noncom_incorrect_mean_gather_2nd_normalised = normalise_custom(response_time_noncom_incorrect_mean_gather_2nd, max(all_responses),min(all_responses));
response_time_noncom_incorrect_std_gather_2nd_normalised = normalise_std(response_time_noncom_incorrect_std_gather_2nd, max(all_responses),min(all_responses));

response_time_noncom_correct_mean_gather_2nd_incorrect_norm = normalise_custom(response_time_noncom_correct_mean_gather_2nd_incorrect, max(all_responses),min(all_responses));
response_time_noncom_correct_std_gather_2nd_incorrect_norm = normalise_std(response_time_noncom_correct_std_gather_2nd_incorrect, max(all_responses),min(all_responses));

response_time_noncom_incorrect_mean_gather_2nd_incorrect_norm = normalise_custom(response_time_noncom_incorrect_mean_gather_2nd_incorrect, max(all_responses),min(all_responses));
response_time_noncom_incorrect_std_gather_2nd_incorrect_norm = normalise_std(response_time_noncom_incorrect_std_gather_2nd_incorrect, max(all_responses),min(all_responses));


figure2=figure;
axes1 = axes('Parent',figure2);
errorbar(response_time_noncom_correct_mean_gather_2nd_normalised,...
    response_time_noncom_correct_std_gather_2nd_normalised...
    ,'LineWidth', 4,'DisplayName',...
    '1^{st} correct, 2^{nd} correct','Color',[0.0977    0.0977    0.4375]);
hold on;
errorbar(response_time_noncom_incorrect_mean_gather_2nd_normalised,...
    response_time_noncom_incorrect_std_gather_2nd_normalised...
    ,'LineWidth', 4,'DisplayName',...
    '1^{st} correct, 2^{nd} error','Color',[0.6953    0.1328    0.1328]);
errorbar(response_time_noncom_correct_mean_gather_2nd_incorrect_norm,...
    response_time_noncom_correct_std_gather_2nd_incorrect_norm...
    ,'LineWidth', 4,'DisplayName',...
    '1^{st} error, 2^{nd} correct','Color',[0.0977    0.0977    0.4375],'LineStyle','-.');

errorbar(response_time_noncom_incorrect_mean_gather_2nd_incorrect_norm,...
    response_time_noncom_incorrect_std_gather_2nd_incorrect_norm...
    ,'LineWidth', 4,'DisplayName',...
    '1^{st}  error, 2^{nd} error','Color',[0.6953    0.1328    0.1328],'LineStyle','-.');

xlim(axes1,[1 6]);
ylim([0 1])
set(axes1,'FontSize',20,'XTickLabel',...
    {'0','3.2','6.4','12.8','25.6','51.2'});

xlabel(['Evidence quality ',char(949), ' (%)']);
ylabel('Normalized response time of 2^{nd} decision');
pubgraph(figure2,28,4,'w')
title_string = 'Response times of 2^{nd} decision detailed';
if(legends)
    legend('show');
end
if(titles)
    title(title_string);
end

set(gcf, 'Position', [161 291 640 680])

if(export)
export_path = [figures_path  'kiani_detailed'];
export_fig (export_path, '-nofontswap', '-linecaps','-png', '-transparent','-r300','-q101', '-cmyk','-painters');
savefig(export_path);
end


end
function coherences = get_coherence_levels(dynamics_and_results)


coherences = zeros(1,size(dynamics_and_results,1));
for i=1:size(dynamics_and_results)
    coherences(i) = dynamics_and_results(i).coherence_level;
end

coherences = unique(coherences);

return
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
                is_motor_correct && dynamics_and_results(i).is_motor_com...
                && dynamics_and_results(i).motor_decision_made == ...
                is_motor_com)
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


function [normalised_vector] = normalise_custom_colour(vector,max,min)
vector = vector';
normalised_vector = zeros(1,size(vector,2));

for i=1:size(vector,2)
    normalised_vector(i) = (vector(i)-min)/(max-min);
end
normalised_vector = normalised_vector';
return;
end
