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

% close all
global legends;
global titles;
global export;
global figures_path;

figures_path = './figures/Fixed_duration';

legends = true;
titles = false;
export = false;


is_first_correct= true;
normalised = true;



stim_onset = dynamics_and_results(1).stim_onset/dynamics_and_results(1).timestep_size;

integration_onset = (dynamics_and_results(1).stim_onset/dynamics_and_results(1).timestep_size) +500;


any = false;
is_motor_com = true;
is_motor_correct = true;
coherence=3.2;
is_late_com = false;
[s_1,s_2,t]=plot_random_trial(dynamics_and_results, coherence, any, is_motor_com, is_motor_correct,is_late_com);


plot_nullcline(x_null_x,x_null_y,y_null_x,y_null_y,s_1,s_2,stim_onset,integration_onset-50)

plot_nullcline(xnull_mc_x,xnull_mc_y,ynull_mc_x,ynull_mc_y,s_1,s_2,integration_onset,t+180)
 
plot_nullcline(x_null_x,x_null_y,y_null_x,y_null_y,s_1,s_2,t+180,8000)


function plot_nullcline(xnull_x,xnull_y,ynull_x,ynull_y,s_1,s_2,start,ending)

figure;
plot(xnull_x,xnull_y)
hold on;
plot(ynull_x,ynull_y)

xlabel('S_1')
ylabel('S_2')


plot(s_1(start:50:ending),s_2(start:50:ending),'.','MarkerSize',15,'Color',[0 0 0]);

v_temp = [0,0;1,1];
plot(v_temp(:,1), v_temp(:,2), 'k--','Color',[0.5 0.5 0.5])

xlim([0 1])
ylim([0 1])
pubgraph(gcf,28,4,'w')
set(gcf, 'Position', [457 311 560 420])

end
%plot_phase_plane(s1,s2,y_1,y_2)

function [s_1,s_2,t]= plot_random_trial(dynamics_and_results, coherence, any, is_motor_com, is_motor_correct,is_late_com)

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
                && dynamics_and_results(i).is_motor_com == is_motor_com && dynamics_and_results(i).is_late_com ==is_late_com)
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

y_5_averaged = calculate_time_window(y_5,100,10);
y_6_averaged = calculate_time_window(y_6,100,10);

plot_time_window(dynamics_and_results, index_to_plot, y_1_averaged,y_2_averaged,100,10,y_5_averaged,y_6_averaged, is_late_com);
[~, t]= max(y_mc_hu);
[s_1,s_2] = convert_y_to_s(y_1,y_2);

return
end

function plot_time_window(dynamics_and_results,index_to_plot,activity_vector1,activity_vector2,time_wind,slide_wind,y_5,y_6,is_late_com)
global legends;
global titles;
global export;
global figures_path;

%this shouldn't be static
dt = dynamics_and_results(1).timestep_size;
b = find(abs(43*(y_5-y_6))>=750);

temp = b(1)+10;

onset = 280;


temp1 = temp +10;
onset = onset+10;
figure;
subplot(3,1,1)
plot((dt*time_wind:dt*slide_wind:dt*(((temp1*slide_wind)+9)-time_wind)),activity_vector1(1:temp-9));hold on;
plot((dt*time_wind:dt*slide_wind:dt*(((temp1*slide_wind)+9)-time_wind)),activity_vector2(1:temp-9));hold on;
% line1=get(gca,'ylim');
% plot([stim_onset stim_onset],line1)
% plot([stim_offset stim_offset],line1)

xlim([(((onset*slide_wind)+9)-time_wind)/2 (((temp1*slide_wind)+9)-time_wind)/2])
ylabel('Firing rate (Hz)');
box off
if(legends)
    legend('Population 1','Population 2');
end
if(titles)
    title('Timecourse of single-trial firing rates (using time window)');
end

% xticklabels({'-3\pi','-2\pi','-\pi','0','\pi','2\pi','3\pi'})
% yticks([-1 -0.8 -0.2 0 0.2 0.8 1])

xticks([])


subplot(3,1,2)
plot((dt*time_wind:dt*slide_wind:dt*(((temp1*slide_wind)+9)-time_wind)),y_5(1:temp-9));hold on;
plot((dt*time_wind:dt*slide_wind:dt*(((temp1*slide_wind)+9)-time_wind)),y_6(1:temp-9));hold on;
% line1=get(gca,'ylim');
% plot([stim_onset stim_onset],line1)
xlim([(((onset*slide_wind)+9)-time_wind)/2  (((temp1*slide_wind)+9)-time_wind)/2])
box off
xticks([])
ylabel('Firing rate (Hz)');

subplot(3,1,3)
plot((dt*time_wind:dt*slide_wind:dt*(((temp1*slide_wind)+9)-time_wind)),43*(y_6(1:temp-9) - y_5(1:temp-9)),'Color',[0 0 0]);hold on;
% line1=get(gca,'ylim');
% plot([stim_onset stim_onset],line1)
xlim([(((onset*slide_wind)+9)-time_wind)/2  (((temp1*slide_wind)+9)-time_wind)/2])
xticks([1406 1724])
xticklabels([0, 1724])
xlabel('Time (ms)'); ylabel('X position');

ylim([-100, 750])
yticks([0, 750])
pubgraph(gcf,28,4,'w')


if(titles)
    title('Difference of firing rates in a single-trial (using time window)');
end


set(gcf, 'Position', [457 311 560 840])

if(export)
export_path = [figures_path 'trial_sample'];
export_fig (export_path, '-nofontswap', '-linecaps','-png', '-transparent','-r300','-q101', '-cmyk','-painters');
savefig(export_path);
end




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
            dynamics_and_results(i).is_motor_correct == is_motor_correct)
        counter = counter+1;
    end
end

y_mc_hu_gather = zeros(counter,trial_length);
y_mc_lu_gather = zeros(counter,trial_length);

j=1;
for i = 1:size(dynamics_and_results,1)
    if(dynamics_and_results(i).coherence_level == coherence &&...
            dynamics_and_results(i).is_motor_com == is_motor_com && ...
            dynamics_and_results(i).is_motor_correct == is_motor_correct)
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
function coherences = get_coherence_levels(dynamics_and_results)


coherences = zeros(1,size(dynamics_and_results,1));
for i=1:size(dynamics_and_results)
    coherences(i) = dynamics_and_results(i).coherence_level;
end

coherences = unique(coherences);

return
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
