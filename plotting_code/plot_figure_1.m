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

figures_path = '../figures_output/Figure1/';




legends = true;
titles = false;
export = false;

is_first_correct= true;
normalised = true;


is_motor_correct = true;
is_motor_com = false;

[y_1, y_2, y_mc_hu, y_mc_lu, y_5, y_6] = ...
    calculate_mean_activity(dynamics_and_results, 3.2, ...
    is_motor_correct, is_motor_com);

[y_1_25, y_2_25, y_mc_hu_25, y_mc_lu_25, y_5_25, y_6_25] = ...
    calculate_mean_activity(dynamics_and_results, 25.6, ...
    is_motor_correct, is_motor_com);

is_single_trial = false;


plot_main_panel_modified(dynamics_and_results,y_1,y_2,y_5,y_6,y_mc_hu,y_mc_lu,y_1_25,y_2_25,y_5_25,y_6_25,y_mc_hu_25,y_mc_lu_25)

plot_accuracies_and_pcom(dynamics_and_results,is_first_correct,normalised);

function plot_main_panel_modified(dynamics_and_results,y_1,y_2,y_5,y_6,y_mc_hu,y_mc_lu,y_1_25,y_2_25,y_5_25,y_6_25,y_mc_hu_25,y_mc_lu_25)

global legends;
global titles;
global export;
global figures_path;

large_left= [char(949), ' 25.6, ', 'L'];
large_right = [char(949), ' 25.6, ', 'R'];

small_left= [char(949), ' 3.2, ', 'L'];
small_right = [char(949), ' 3.2, ', 'R'];


large = [char(949), ' 25.6',];
small = [char(949),' 3.2',];


response_time_for_25 = find(y_6_25>17.4,1);

response_time_for_32 = find(y_6>17.4,1);



figure1=figure;
subplot(3,1,1);
title('Sensorimotor');
hold on;
h = zeros(1,6);
h(1) = plot((1400:2:response_time_for_25), y_1_25(1400:2:response_time_for_25),'Color',[0.1172    0.5625    1.0000],'DisplayName',large_left);
h(2) = plot((1400:2:response_time_for_25), y_2_25(1400:2:response_time_for_25),'Color',[1.0000    0.4961    0.3125],'DisplayName',large_right);

%%% grey stuff
h(3) = plot((response_time_for_25:2:5000), y_1_25(response_time_for_25:2:5000),'Color',[0.5000    0.5000    0.5000],'DisplayName',large_left);
h(4) = plot((response_time_for_25:2:5000), y_2_25(response_time_for_25:2:5000),'Color',[0.5000    0.5000    0.5000],'DisplayName',large_right);


h(5) = plot((1400:2:response_time_for_32), y_1(1400:2:response_time_for_32),'LineStyle','-.','Color',[0.1172    0.5625    1.0000],'DisplayName',small_left);
h(6) = plot((1400:2:response_time_for_32), y_2(1400:2:response_time_for_32),'LineStyle','-.','Color',[1.0000    0.4961    0.3125],'DisplayName',small_right);

%%% grey stuff
h(7) = plot((response_time_for_32:2:5000),y_1(response_time_for_32:2:5000),'Color',[0.5000    0.5000    0.5000],'DisplayName',large_left);
h(8) = plot((response_time_for_32:2:5000),y_2(response_time_for_32:2:5000),'Color',[0.5000    0.5000    0.5000],'DisplayName',large_right);

xlim([1500 5000])
xticks([])
xticklabels({})
if(legends)
    legend([h(1),h(2),h(5),h(6)]);
end
if(titles)
   title(['Decision' title_string]);
end

yticks([0, 30])
yticklabels({'0','40'})
ylabel('Firing rate (Hz)');
%xlabel('Time (ms)');

subplot(3,1,2);
title('Uncertainty Monitoring');
hold on;

v = zeros(1,6);
v(1) = plot((1500:2:response_time_for_25),y_mc_hu_25(1500:2:response_time_for_25),'Color',[1 0 0],'DisplayName',large);
v(2) = plot((1500:2:response_time_for_32), y_mc_hu(1500:2:response_time_for_32),'LineStyle','-.','Color',[1 0 0],'DisplayName',small);

%%% grey stuff
v(3) = plot((response_time_for_25:2:5000), y_mc_hu_25(response_time_for_25:2:5000),'Color',[0.5000    0.5000    0.5000],'DisplayName',large_left);
v(4) = plot((response_time_for_32:2:5000), y_mc_hu(response_time_for_32:2:5000),'Color',[0.5000    0.5000    0.5000],'DisplayName',large_right);


v(5) = plot((1500:2:response_time_for_32), y_mc_lu(1500:2:response_time_for_32),'LineStyle','-.','Color',[0.4    0.73        0.27],'DisplayName',small);
v(6) = plot((1500:2:response_time_for_25), y_mc_lu_25(1500:2:response_time_for_25),'Color',[0.4    0.73        0.27],'DisplayName',large);

%%% grey stuff
v(7) = plot((response_time_for_25:2:5000), y_mc_lu_25(response_time_for_25:2:5000),'Color',[0.5000    0.5000    0.5000],'DisplayName',large_left);
v(8) = plot((response_time_for_32:2:5000), y_mc_lu(response_time_for_32:2:5000),'Color',[0.5000    0.5000    0.5000],'DisplayName',large_right);


xlim([1500,5000]);
xticks([])
xticklabels({})
if(legends)
    legend([v(1),v(2),v(5),v(6)]);
end
if(titles)
    title(['Uncertainty monitoring' title_string]);
end

ylabel('Firing rate (Hz)');
%xlabel('Time (ms)');


subplot(3,1,3);
title('Motor');
hold on;
g = zeros(1,8);
g(1) = plot((1500:2:response_time_for_25), y_5_25(1500:2:response_time_for_25),'Color',[0.1172    0.5625    1.0000]);
g(2) = plot((1500:2:response_time_for_25), y_6_25(1500:2:response_time_for_25),'Color',[1.0000    0.4961    0.3125]);

%%% grey stuff
g(3) = plot((response_time_for_25:2:5000), y_5_25(response_time_for_25:2:5000),'Color',[0.5000    0.5000    0.5000],'DisplayName',large_left);
g(4) = plot((response_time_for_25:2:5000), y_6_25(response_time_for_25:2:5000),'Color',[0.5000    0.5000    0.5000],'DisplayName',large_right);

g(5)= plot((1500:2:response_time_for_32), y_5(1500:2:response_time_for_32),'LineStyle','-.','Color',[0.1172    0.5625    1.0000]);
g(6)= plot((1500:2:response_time_for_32), y_6(1500:2:response_time_for_32),'LineStyle','-.','Color',[1.0000    0.4961    0.3125]);

%%% grey stuff
g(7)= plot((response_time_for_32:2:5000), y_5(response_time_for_32:2:5000),'Color',[0.5000    0.5000    0.5000],'DisplayName',large_left);
g(8)= plot((response_time_for_32:2:5000), y_6(response_time_for_32:2:5000),'Color',[0.5000    0.5000    0.5000],'DisplayName',large_right);


ylabel('Firing rate (Hz)');
xlabel('Time (ms)');
xlim([1500,5000]);
xticks([1500, 5000])
xticklabels({'0', '2500'})

if(legends)
    %legend('Left','Right');
end
if(titles)
    title(['Motor output' title_string]);
end
pubgraph(figure1,28,4,'w');

if(export)
export_path = [figures_path  'new_panel'];
export_fig (export_path, '-nofontswap', '-linecaps','-png', '-transparent','-r300','-q101', '-cmyk','-painters');
savefig(export_path);
end


%%% Threshold line:
x=1500:2:5000;
y=17.4;
plot(x,y*ones(size(x)),'LineStyle','--','Color',[0 0 0]);


set(gcf, 'Position', [175 71 560 900])

end


function[y_1_mean, y_2_mean, y_mc_hu_mean, y_mc_lu_mean, y_5_mean,...
    y_6_mean] = ...
    calculate_mean_activity(dynamics_and_results, coherence, ...
    is_motor_correct, is_motor_com)

coherences = get_coherence_levels(dynamics_and_results);

if(~ismember(coherence,coherences))
    error('coherence level not found in data');
end

counter = 0;
trial_length =  dynamics_and_results(1).trial_length;

%get number of each type, before initialising vectors.


for i = 1:size(dynamics_and_results,1)
    if(dynamics_and_results(i).coherence_level == coherence && dynamics_and_results(i).is_motor_correct == is_motor_correct && dynamics_and_results(i).is_motor_com == is_motor_com)
        counter = counter+1;
    end
end

y_1_gather = zeros(counter,trial_length);
y_2_gather = zeros(counter,trial_length);
y_mc_hu_gather = zeros(counter,trial_length);
y_mc_lu_gather = zeros(counter,trial_length);
y_5_gather = zeros(counter,trial_length);
y_6_gather = zeros(counter,trial_length);

j=1;
for i = 1:size(dynamics_and_results,1)
    if(dynamics_and_results(i).coherence_level == coherence && dynamics_and_results(i).is_motor_correct == is_motor_correct && dynamics_and_results(i).is_motor_com == is_motor_com)
        y_1_gather(j,:)=dynamics_and_results(i).y_1;
        y_2_gather(j,:)=dynamics_and_results(i).y_2;
        y_mc_hu_gather(j,:)=dynamics_and_results(i).y_mc_hu;
        y_mc_lu_gather(j,:)=dynamics_and_results(i).y_mc_lu;
        y_5_gather(j,:)=dynamics_and_results(i).y_5;
        y_6_gather(j,:)=dynamics_and_results(i).y_6;
        j=j+1;
    end
end


% Non CoM Correct Mean
y_1_mean(1,:) = mean(y_1_gather);
y_2_mean(1,:) = mean(y_2_gather);
y_mc_hu_mean(1,:) = mean(y_mc_hu_gather);
y_mc_lu_mean(1,:) = mean(y_mc_lu_gather);
y_5_mean(1,:) = mean(y_5_gather);
y_6_mean(1,:) = mean(y_6_gather);

return;
end



function coherences = get_coherence_levels(dynamics_and_results)


coherences = zeros(1,size(dynamics_and_results,1));
for i=1:size(dynamics_and_results)
    coherences(i) = dynamics_and_results(i).coherence_level;
end

coherences = unique(coherences);

return
end


function plot_accuracies_and_pcom(dynamics_and_results, is_first_correct,normalised)
global legends;
global titles;
global export;
global figures_path;
coherences = get_coherence_levels(dynamics_and_results);
n_of_coherences= size(coherences,2);

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
   
    p_correct(i) = (noncom_correct_counter+com_correct_counter)/all_trials;
    p_correct_noncom(i) = (noncom_correct_counter)/all_trials;
    p_nondecision(i)= nondecision_counter/all_trials;
    
    [initiation_time_noncom_correct_mean_gather(i),initiation_time_noncom_correct_std_gather(i),...
        response_time_noncom_correct_mean_gather(i),response_time_noncom_correct_std_gather(i)] = ...
        calculate_mean_times(dynamics_and_results, coherence, ...
        true, false, is_first_correct, coupled);
    
    [initiation_time_noncom_incorrect_mean_gather(i),initiation_time_noncom_incorrect_std_gather(i),...
        response_time_noncom_incorrect_mean_gather(i),response_time_noncom_incorrect_std_gather(i)] = ...
        calculate_mean_times(dynamics_and_results, coherence, ...
        false, false, is_first_correct, coupled);
    
    [initiation_time_com_correct_mean_gather(i),initiation_time_com_correct_std_gather(i),...
        response_time_com_correct_mean_gather(i),response_time_com_correct_std_gather(i)] = ...
        calculate_mean_times(dynamics_and_results, coherence, ...
        true, true, is_first_correct, coupled);
    
    [initiation_time_com_incorrect_mean_gather(i),initiation_time_com_incorrect_std_gather(i),...
        response_time_com_incorrect_mean_gather(i),response_time_com_incorrect_std_gather(i)] = ...
        calculate_mean_times(dynamics_and_results, coherence, ...
        false, true, is_first_correct, coupled);
 
    
end


%%%%%%%%Fig1: Choice accuracy


% Set up fittype and options.
ft = fittype( '1-0.5*exp(-(x/a).^b)', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [0.68988590882247 0.217721776227971];

% Fit model to data.
[fitresult, gof] = fit( [0,3.2,6.4,12.8,25.6,51.2]', p_correct', ft, opts );

%%% extract data from fitted curve
figuretemp = figure;
plot(fitresult);
h = findobj(figuretemp,'Type','line');
x=get(h,'Xdata');
y=get(h,'Ydata');
close;
%%%%%%%


figure3=figure;
plot(x,y,'LineWidth', 4,'Color',[0 0 0]);
hold on;
scatter([0,3.2,6.4,12.8,25.6,51.2],p_correct,'o','LineWidth',4,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor', [0 0 0]);
xlabel(['Evidence quality ',char(949), ' (%)']);
ylabel('Probability correct');
xlim([0,51.2])
ylim([min(y) 1])
pubgraph(figure3,28,4,'w')
title_string = 'Choice accuracy';


if(legends)
    %legend('show');
end
if(titles)
    title(title_string);
end

if(export)
export_path = [figures_path  'choice_accuracy'];
export_fig (export_path, '-nofontswap', '-linecaps','-png', '-transparent','-r300','-q101', '-cmyk','-painters');
savefig(export_path);
end

%%%%%%%%

%%%%%%%%Figure2: initiation time

figure4=figure;
axes1 = axes('Parent',figure4);
hold on;


normalised_initiation_time_noncom_correct = normalise_custom(initiation_time_noncom_correct_mean_gather, max(initiation_time_noncom_incorrect_mean_gather),min(initiation_time_noncom_correct_mean_gather));
normalised_initiation_time_noncom_correct_std = normalise_std(initiation_time_noncom_correct_std_gather, max(initiation_time_noncom_incorrect_mean_gather),min(initiation_time_noncom_correct_mean_gather));



normalised_initiation_time_noncom_incorrect = normalise_custom(initiation_time_noncom_incorrect_mean_gather, max(initiation_time_noncom_incorrect_mean_gather),min(initiation_time_noncom_correct_mean_gather));
normalised_initiation_time_noncom_incorrect_std = normalise_std(initiation_time_noncom_incorrect_std_gather, max(initiation_time_noncom_incorrect_mean_gather),min(initiation_time_noncom_correct_mean_gather));

errorbar(normalised_initiation_time_noncom_correct,...
    normalised_initiation_time_noncom_correct_std...
    ,'LineWidth', 4,'DisplayName',...
    'correct','Color',[0 0 1]);
errorbar(normalised_initiation_time_noncom_incorrect,...
    normalised_initiation_time_noncom_incorrect_std...
    ,'LineWidth', 4,'DisplayName',...
    'error','Color',[1 0 0]);

% scatter([0,3.2,6.4,12.8,25.6,51.2],normalised_initiation_time_noncom_correct,'o','LineWidth',4,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor', [0 0 1]);
% scatter([0,3.2,6.4,12.8,25.6,51.2],normalised_initiation_time_noncom_incorrect,'o','LineWidth',4,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor', [1 0 0]);

xlim(axes1,[1 6]);
ylim([0 1])
set(axes1,'FontSize',20,'XTickLabel',...
    {'0','3.2','6.4','12.8','25.6','51.2'});

xlabel(['Evidence quality ',char(949), ' (%)']);
ylabel('Normalized response time');
pubgraph(figure4,28,4,'w')
title_string = 'Initiation times v.s coherence level';
if(legends)
    legend('show');
end
if(titles)
    title(title_string);
end

if(export)
export_path = [figures_path  'reaction_times'];
export_fig (export_path, '-nofontswap', '-linecaps','-png', '-transparent','-r300','-q101', '-cmyk','-painters');
savefig(export_path);
end

%%%%%%%%

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