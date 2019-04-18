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

% clearvars -except dynamics_and_results dynamics_and_results_2nddecision;

global legends;
global titles;
global export;
global figures_path;

figures_path = '../figures_output/supp';

legends = true;
titles = false;
export = false;


is_first_correct=    true;
normalised = true;

any = false;
is_motor_com = false;
is_motor_correct = false;
coherence=3.2;
is_late_com = false;
onresponse=true;
plot_random_trial(dynamics_and_results, coherence, any, is_motor_com, is_motor_correct,is_late_com,onresponse)

function plot_random_trial(dynamics_and_results, coherence, any, is_motor_com, is_motor_correct,is_late_com,onresponse)

coherences = get_coherence_levels(dynamics_and_results);

if(~ismember(coherence,coherences))
    error('coherence level not found in data');
end

counter = 0;
%get number of each type, before initialising vectors.


for i = 1:size(dynamics_and_results,1)
    if(dynamics_and_results(i).coherence_level == coherence && dynamics_and_results(i).is_motor_correct == is_motor_correct && dynamics_and_results(i).is_motor_com == is_motor_com && dynamics_and_results(i).is_late_com == is_late_com)
        if(is_motor_com ==  dynamics_and_results(i).is_motor_com && is_motor_com == check_com_advanced_condition(dynamics_and_results, i))
            counter = counter+1;
        end
    end
end

indecies = zeros(1,counter);
j=1;

if(any)
    for i = 1:size(dynamics_and_results,1)
        if(dynamics_and_results(i).coherence_level == coherence)
            indecies(j) = i;
            j=j+1;
        end
    end
else
    for i = 1:size(dynamics_and_results,1)
        if(dynamics_and_results(i).coherence_level == coherence &&...
                dynamics_and_results(i).is_motor_correct == is_motor_correct...
                && dynamics_and_results(i).is_motor_com == ...
                is_motor_com &&...
                dynamics_and_results(i).is_late_com ==is_late_com &&...
                dynamics_and_results(i).motor_decision_made)
            if(is_motor_com ==  dynamics_and_results(i).is_motor_com && is_motor_com == check_com_advanced_condition(dynamics_and_results, i))
                indecies(j) = i;
                j=j+1;
            end
        end
    end
end

%pick a random element

index_to_plot= indecies(randi(numel(indecies)));

y_1 = dynamics_and_results(index_to_plot).y_1;
y_2 = dynamics_and_results(index_to_plot).y_2;
y_mc_hu = dynamics_and_results(index_to_plot).y_mc_hu;
y_mc_lu = dynamics_and_results(index_to_plot).y_mc_lu;
y_5 = dynamics_and_results(index_to_plot).y_5;
y_6 = dynamics_and_results(index_to_plot).y_6;



y_1_averaged = calculate_time_window(y_1,100,10);
y_2_averaged = calculate_time_window(y_2,100,10);
y_mc_hu_averaged = calculate_time_window(y_mc_hu,100,10);
y_5_averaged = calculate_time_window(y_5,100,10);
y_6_averaged = calculate_time_window(y_6,100,10);

% plot_main_panel(y_1,y_2,y_5,y_6,y_mc_hu,y_mc_lu,true,coherence);

plot_time_window(dynamics_and_results, index_to_plot, y_1_averaged,y_2_averaged,y_mc_hu_averaged,y_5_averaged,y_6_averaged,100,10);

end

function plot_time_window(dynamics_and_results,index_to_plot,activity_vector1,activity_vector2,y_mc_hu_averaged,y_5_averaged,y_6_averaged,time_wind,slide_wind)
global legends;
global titles;
global export;
global figures_path;

%this shouldn't be static
dt = dynamics_and_results(index_to_plot).timestep_size;
trial_length = dynamics_and_results(index_to_plot).trial_length;

response = dynamics_and_results(index_to_plot).response_time;
movement = dynamics_and_results(index_to_plot).response_time + dynamics_and_results(index_to_plot).movement_duration;
stim_onset = dynamics_and_results(index_to_plot).stim_onset;

averaging_t = (dt*time_wind:dt*slide_wind:dt*(trial_length-time_wind));

movement_avg = find(averaging_t>=movement,1);

figure;
subplot(3,1,1)
plot((dt*time_wind:dt*slide_wind:dt*(trial_length-time_wind)),activity_vector1(1:end-10),'Color',[0.1172    0.5625    1.0000]);hold on;
plot((dt*time_wind:dt*slide_wind:dt*(trial_length-time_wind)),activity_vector2(1:end-10),'Color',[1.0000    0.4961    0.3125]);hold on;
ylabel('Firing rate (Hz)');
title('Sensorimotor')
box off
% xlim([stim_onset movement])
xticks([])
if(legends)
    legend('Left','Right');
end
if(titles)
    title('Timecourse of single-trial firing rates (using time window)');
end

subplot(3,1,2)
plot((dt*time_wind:dt*slide_wind:dt*(trial_length-time_wind)),y_mc_hu_averaged(1:end-10),'Color',[1 0 0]);hold on;
ylabel('Firing rate (Hz)');
title('Uncertainty Monitoring');
box off
% xlim([stim_onset movement])
xticks([])


subplot(3,1,3)
plot((dt*time_wind:dt*slide_wind:dt*(trial_length-time_wind)),y_5_averaged(1:end-10),'Color',[0.1172    0.5625    1.0000]);hold on;
plot((dt*time_wind:dt*slide_wind:dt*(trial_length-time_wind)),y_6_averaged(1:end-10),'Color',[1.0000    0.4961    0.3125]);hold on;
xticks([])
yticks()

xticks([stim_onset movement])
xticklabels({'0', num2str(movement)})
yticklabels([0 17])
% xlim([stim_onset movement])
ylabel('Firing rate (Hz)');
title('Motor')
box off
xlabel('Time (ms)'); ylabel('Firing rate (Hz)');
pubgraph(gcf,28,4,'w');

set(gcf, 'Position', [303 76 581 895])
if(titles)
    title('Difference of firing rates in a single-trial (using time window)');
end

set(gcf, 'Position', [457 311 560 840])

if(export)
export_path = [figures_path + 'time_window'];
export_fig (export_path, '-nofontswap', '-linecaps','-png', '-transparent','-m10','-q101', '-cmyk','-painters');
savefig(export_path);
end


end


function[time_window_averaged] = ...
    calculate_time_window(activity_vector, time_wind, slide_wind)


trial_length = size(activity_vector,2);

time_window_averaged = zeros(1,(trial_length/slide_wind)-9);
time_window_averaged(1) = (mean(activity_vector(1:time_wind))) ;



for t = 2:(trial_length-time_wind)/slide_wind
    time_window_averaged(t) = ...
        (mean(activity_vector(slide_wind*t:slide_wind*t+time_wind))) ;
end


end

%% helper functions:

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

function [itiscom] = check_com_advanced_condition(dynamics_and_results, index)

itiscom = false;

x_traj = dynamics_and_results(index).y_6-dynamics_and_results(index).y_5;


if(dynamics_and_results(index).is_motor_correct)
    x_traj = dynamics_and_results(index).y_5-dynamics_and_results(index).y_6;
end

smoothed_trajectory = filter(ones(1,50)/50,1,x_traj);
delta_of_trajectory = diff(sign(smoothed_trajectory));

points_of_change = find(delta_of_trajectory);

if(size(points_of_change,2)<2)
    return;
end

for i=size(points_of_change,2):-1:2
    if(max(smoothed_trajectory(points_of_change(i-1):points_of_change(i)))>0)
        sign_of_response = sign(smoothed_trajectory((dynamics_and_results(index).movement_duration + dynamics_and_results(index).response_time)/dynamics_and_results(index).timestep_size));
        sign_of_max_point_in_bump =  sign(max(smoothed_trajectory(points_of_change(i-1):points_of_change(i))));
        
        
        if(sign_of_response *  sign_of_max_point_in_bump == -1)
            if(max(abs(smoothed_trajectory(points_of_change(i-1):points_of_change(i))))>2)
                itiscom=true;
                return;
            end

        end

    end
    
end

return

end
