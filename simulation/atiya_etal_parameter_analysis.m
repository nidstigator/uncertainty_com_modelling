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


global dt; global trial_length; global stim_onset; global stim_duration;

global j_mc0;

global mu;


global inhibitory_onset;
global uncertainty_onset;



dt = 0.5;
trial_length = 4000/dt;
stim_onset = 900;
stim_duration = 850;

%%param analysis variables - default values!
mu = 30;
j_mc0 = 0.009;

uncertainty_onset=500;
inhibitory_onset=500;


%%%%Threshold values:
decision_module_threshold = 35.5;
motor_target_threshold = 17.4;
%%%%%%%%%%%%%%%%%%%%%%

% prompt = 'Please enter the coherence level, or please enter 100 for all: \n';
coherences = 6.4;

coherence_string = [num2str(coherences) '_coh'];

if(coherences==100)
    coherences = [0,3.2,6.4,12.8,25.6,51.2];
    coherence_string = 'all_coh';
end
clearvars -except coherences inhibitory_onset uncertainty_onset mu j_mc0  ...
    motor_target_threshold decision_module_threshold ...
    stim_duration stim_onset trial_length dt 
close all;
prompt = 'Please enter the number of trials: \n';
n_of_trials = input(prompt);


prompt = '0 for offset testing, 1 for j_mc and mu testing: \n';
testing_mode = input(prompt);

tStart=tic;

dynamics_and_results = repmat(struct('block_no',0,'block_size',0,...
    'trial_no',0,'overall_mc_input',[],'coherence_level',0,'timestep_size',0,'trial_length',0,...
    'stim_onset',0,'stim_duration',0,...
    'y_1',[],'y_2',[],'y_5',[],'y_6',[],...
    'y_mc_hu',[], 'y_mc_lu',[],'c',0, 'decision_module_threshold',0,...
    'motor_target_threshold',0,'motor_maximum_value_left',0,...
    'motor_maximum_value_right',0,'decision_module_crossed',false,...
    'is_motor_correct',false,...
    'response_time',0,'movement_duration',0,'is_motor_com',false,...
    'is_late_com',false,'motor_decision_made',false)...
    , n_of_trials*size(coherences,2), 1);




% coherences for loop
index = 0;
analysis_index=0;

no_of_blocks = size(coherences,2);

number_of_simulations = 3;

parameter_values_to_test_for_j_mc = [0, 0.0005, 0.001, 0.002, 0.005, 0.009, 0.02, 0.05, 0.07];
parameter_values_to_test_for_mu = [20, 23, 26, 29 , 32, 35, 38, 40];

parameter_values_to_test_for_lu_offset = [0,100,200,300,400,500,600];
parameter_values_to_test_for_hu_offset = [0,100,200,300,400,500,600];
param_to_test_1 = parameter_values_to_test_for_j_mc;
param_to_test_2 = parameter_values_to_test_for_mu;
analysis_results = repmat(struct('j_mc_value',0,'constant_input_value',0,'p_correct',0,'p_com',0,'mean_rt',0),(size(param_to_test_1,2)*size(param_to_test_2,2)),1);


if(testing_mode==0)
    param_to_test_1 = parameter_values_to_test_for_lu_offset;
    param_to_test_2 = parameter_values_to_test_for_hu_offset;
    analysis_results = repmat(struct('inhibitory_onset','uncertainty_onset','p_correct',0,'p_com',0,'mean_rt',0),(size(param_to_test_1,2)*size(param_to_test_2,2)),1);
end

for o = 1 :size(param_to_test_1,2)
    
    if(testing_mode==0)
        inhibitory_onset = param_to_test_1(o);
    else
        j_mc0 = param_to_test_1(o);
    end
   
    for q = 1: size(param_to_test_2,2)
        analysis_index = analysis_index + 1;
        
        if(testing_mode==0)
            uncertainty_onset = param_to_test_2(q);
        else
            mu = param_to_test_2(q);
        end
        index = 0;
        for k = 1: no_of_blocks
            coh = coherences(k);
            % Trials for loop.
            for i = 1:n_of_trials
                index = index+1;
                [y_1,y_2,s1,s2,y_mc_hu,y_mc_lu,decision_module_crossed,...
                    response_time,overall_mc_input] = nonlinear_decision_module_and_mc_integrate...
                    (coh,decision_module_threshold);
                

                [y_5,y_6,response_time] = ...
                    hand_module_integrate(y_1,y_2,response_time);
                
                % Set these values in case functions don't execute due to
                % nondecision
                motor_decision_made=false;
                movement_duration = 0;
                is_motor_correct = false;
                is_motor_com = false;
                is_late_com = false;
                %%%%%%%%%%
                
                if(~(isempty(y_5) || isempty(y_6)))
                    [is_motor_correct, motor_decision_made, movement_duration] = ...
                        check_if_motor_target_reached...
                        (y_5,y_6,...
                        motor_target_threshold,response_time);
                    
                    [is_motor_com,is_late_com] = check_com(y_5,y_6,motor_decision_made,motor_target_threshold,is_motor_correct);
                end
                %%%% Save in struct:
                dynamics_and_results(index).block_no=k;
                dynamics_and_results(index).block_size=n_of_trials;
                dynamics_and_results(index).trial_no=i;
                dynamics_and_results(index).coherence_level=coh;
                dynamics_and_results(index).timestep_size=dt;
                dynamics_and_results(index).trial_length=trial_length;
                dynamics_and_results(index).stim_onset=stim_onset;
                dynamics_and_results(index).stim_duration=stim_duration;
                dynamics_and_results(index).y_1=y_1;
                dynamics_and_results(index).y_2=y_2;
                dynamics_and_results(index).y_5=y_5;
                dynamics_and_results(index).y_6=y_6;
                dynamics_and_results(index).y_mc_hu=y_mc_hu;
                dynamics_and_results(index).y_mc_lu=y_mc_lu;
                dynamics_and_results(index).overall_mc_input = overall_mc_input;
                dynamics_and_results(index).decision_module_threshold = decision_module_threshold;
                dynamics_and_results(index).motor_target_threshold = motor_target_threshold;
                dynamics_and_results(index).motor_maximum_value_left = max(y_5);
                dynamics_and_results(index).motor_maximum_value_right = max(y_6);
                dynamics_and_results(index).decision_module_crossed = decision_module_crossed;
                dynamics_and_results(index).motor_decision_made = motor_decision_made;
                dynamics_and_results(index).is_motor_correct = is_motor_correct;
                dynamics_and_results(index).movement_duration = movement_duration * dt;
                dynamics_and_results(index).response_time = response_time * dt;
                dynamics_and_results(index).is_motor_com = is_motor_com;
                dynamics_and_results(index).is_late_com= is_late_com;
                %%%%%%%%%%%
                
                %%%% Clear vars before next trial
                clearvars -except analysis_results dynamics_and_results k i o q dt trial_length stim_onset...
                    stim_duration stim_offset decision_module_threshold ...
                    motor_preparation_module_threshold motor_target_threshold ...
                    hand_module_delay mc_module_delay click_delay is_fixed_duration...
                    coh n_of_trials tStart filename coherences index across_trial_effect no_of_blocks param_to_test_1 param_to_test_2 uncertainty_onset inhibitory_onset j_mc0 mu analysis_index testing_mode
                %%%%
            end
            
        end

        [p_com,p_correct,p_com_correct,p_com_incorrect,p_late_com_correct,p_late_com_incorrect] = calculate_accuracy_and_pcom(dynamics_and_results);
        [mean_rt] = calculate_mean_rt(dynamics_and_results);
        if(testing_mode==0)
            analysis_results(analysis_index).uncertainty_onset = uncertainty_onset;
            analysis_results(analysis_index).inhibitory_onset = inhibitory_onset;
        else
            analysis_results(analysis_index).j_mc_value = j_mc0;
            analysis_results(analysis_index).constant_input_value = mu;
        end
        analysis_results(analysis_index).p_correct =p_correct;
        analysis_results(analysis_index).p_com =p_com;
        analysis_results(analysis_index).mean_rt = mean_rt;
        analysis_results(analysis_index).p_com_correct = p_com_correct;
        analysis_results(analysis_index).p_com_incorrect = p_com_incorrect;
        analysis_results(analysis_index).p_late_com_correct = p_late_com_correct;
        analysis_results(analysis_index).p_late_com_incorrect = p_late_com_incorrect;
        
        
        dynamics_and_results = repmat(struct('block_no',0,'block_size',0,...
            'trial_no',0,'overall_mc_input',[],'coherence_level',0,'timestep_size',0,'trial_length',0,...
            'stim_onset',0,'stim_duration',0,'stim_offset',0,'is_fixed_duration',false,...
            'y_1',[],'y_2',[],'y_5',[],'y_6',[],...
            'y_mc_hu',[], 'y_mc_lu',[],'c',0, 'decision_module_threshold',0,...
            'motor_target_threshold',0,'motor_maximum_value_left',0,...
            'motor_maximum_value_right',0,'decision_module_crossed',false,...
            'is_motor_correct',false,...
            'response_time',0,'movement_duration',0,'is_motor_com',false,...
            'is_late_com',false,'motor_decision_made',false)...
            , n_of_trials*size(coherences,2), 1);

        clearvars  -except dynamics_and_results analysis_results k o q dt trial_length stim_onset...
            stim_duration stim_offset decision_module_threshold ...
            motor_target_threshold ...
            coh n_of_trials tStart filename coherences across_trial_effect no_of_blocks param_to_test_1 param_to_test_2 uncertainty_onset inhibitory_onset j_mc0 mu analysis_index testing_mode
    end
end


clearvars -except analysis_results tStart filename

tElapsed=toc(tStart);

fprintf('Time taken to run this: %f \n',tElapsed)

clearvars tStart tElapsed;


%old save function doesn't work with large variables.
%save(filename,'-v7.3');

function [y_1,y_2,s1,s2,y_mc_uncertainty,y_mc_inhibitory,decision_module_crossed,response_time,overall_mc_input] = nonlinear_decision_module_and_mc_integrate(coh,decision_module_threshold)
global inhibitory_onset; global uncertainty_onset; global j_mc0; global mu;
global dt; global trial_length; global stim_onset;

response_time = 0;
%%%%%%%%%%%% DECISION PARAMETERS
%%%% Synaptic time and other constants


tnmda = 100;    % NMDAr
tampa = 2;      % AMPAr
gamma = 0.641;  % Gamma

% FI curve parameters

a = 270; b = 108; d = 0.1540;  % Parameters for excitatory cells

% Parameters to be varied

mu0       = 30.0;      % External stimulus strength
noise_amp = 0.025;      % Noise amplitude into selective populations


%---- Initial conditions and clearing variables -----------
s1_in=0.1; s2_in=0.1;
I_eta1_in = noise_amp*randn ; I_eta2_in = noise_amp*randn ;


% Intialise and vectorise variables to be used in loops below

s1 = s1_in.*ones(1,trial_length); s2 = s2_in.*ones(1,trial_length);
I_eta1 = I_eta1_in.*ones(1,trial_length);
I_eta2 = I_eta2_in.*ones(1,trial_length);
%%%%%%%%%%%%

%%%%%%%%%%%% UNCERTAINTY MODULE PARAMETERS:

gain_uncertainty = 1;     % Gain of input-output function (cf. Wong & Wang, 2006)
tau_mc_uncertainty = 150;    % Membrane time constant (cf. Wilson and Cowan, 1972)
gain_inh = 1;     % Gain of input-output function (cf. Wong & Wang, 2006)
tau_mc_inhibitory = 150;    % Membrane time constant (cf. Wilson and Cowan, 1972)


j_mc_v_inh = 0.5;% Inhibition from inhibitory to uncertainty neural population
j_self_uncertainty = 0;      % Self excitation constant for uncertainty.
%%%%%%%%%%%%

% Initialise state vectors with zeros (for better performance)
y_1 = zeros(1,trial_length);
y_2 = zeros(1,trial_length);
% Initialise state vectors with zeros (for better performance)
y_mc_inhibitory = zeros(1,trial_length);
y_mc_uncertainty = zeros(1,trial_length);
overall_mc_input = zeros(1,trial_length);

Isyn1 = zeros(1,trial_length);
Isyn2 = zeros(1,trial_length);

motor_preparation = 0;
decision_module_crossed = false;
for t=1:trial_length-1
    %Constant effective external current input
    I0E1 = 0.3255; I0E2 = 0.3255;
    
    % External stimulus
    JAext = 0.00052; % Synaptic coupling constant to external inputs
    
    I_stim_1 = (stim_onset/dt<t & ...
        motor_preparation==0)*...
        (JAext * mu0 * (1-coh/100)); % To population 1
    I_stim_2 = (stim_onset/dt<t & ...
        motor_preparation==0)*...
        (JAext * mu0 * (1+coh/100)); % To population 2
    
    
    % Recurrent synaptic coupling constants
    JN11 = 0.2440182353; JN22 = 0.2440182353;
    JN12 = 0.0497; JN21 = 0.0497;
    
    % Resonse function of competiting excitatory population 1
    Isyn1(t) = JN11 .* s1(t) - JN12 .* s2(t) + I0E1 + I_stim_1 ...
        + I_eta1(t)+ (j_mc0 * y_mc_uncertainty(t));

    y_1(t)  = (a .* Isyn1(t) - b) ./ (1 ...
        - exp(-d .* (a .* Isyn1(t)-b)));
    
    % Response function of competiting excitatory population 2
    Isyn2(t) = JN22 .* s2(t) - JN21 .* s1(t) + I0E2 + I_stim_2 ...
        + I_eta2(t) + (j_mc0 * y_mc_uncertainty(t));
    
    y_2(t)  = (a .* Isyn2(t) - b) ./ (1 ...
        - exp(-d.*(a.*Isyn2(t)-b)));
    
    % Dynamical equations
    
    % Mean NMDA-receptor dynamics
    s1(t+1) = s1(t) + dt * (-(s1(t)/tnmda) ...
        + (1 - s1(t)) * gamma * y_1(t)/1000);
    s2(t+1) = s2(t) + dt*(-(s2(t)/tnmda) ...
        + (1 - s2(t)) * gamma * y_2(t)/1000);
    
    % Noise through synaptic currents of pop1 and 2
    I_eta1(t+1) = I_eta1(t) + (dt/tampa) * (-I_eta1(t))...
        + sqrt(dt/tampa) * noise_amp*randn;
    I_eta2(t+1) = I_eta2(t) + (dt/tampa) * (-I_eta2(t))...
        + sqrt(dt/tampa) * noise_amp*randn;
    
    % To ensure firing rates are always positive
    if (y_1(t) < 0 )
        y_1(t) = 0;
    end
    if (y_2(t) < 0)
        y_2(t) = 0;
    end
    
    %to record moment of threshold crossing for motor preparation
    %this is for 'more natural threshold crossing'
    if ((y_1(t) > decision_module_threshold ...
            || y_2(t) > decision_module_threshold) ...
            && ~decision_module_crossed)
        response_time = t;
        motor_preparation =1;
        decision_module_crossed = true;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%UNCERTAINTY MODULE INTEGRATE:
    % Conditional inhibitory control.
    inh_top_down_inhibition = ((stim_onset+inhibitory_onset)/dt>t) * 1000;
    uncertainty_top_down_inhibition = ((stim_onset+uncertainty_onset)/dt>t) * 1000;
    
    custom_top_down_inhibition = (decision_module_crossed) * 3000;
    
    % Total input coming in to uncertainty population (inhibitory and uncertainty)
    i_uncertainty_total = gain_uncertainty * (mu - ...
        j_mc_v_inh * y_mc_inhibitory(t)) - inh_top_down_inhibition - custom_top_down_inhibition;
    
    i_inhibitory_total = gain_inh * (j_mc0 * (y_1(t) +...
        y_2(t)) + j_self_uncertainty * y_mc_inhibitory(t)) - uncertainty_top_down_inhibition - custom_top_down_inhibition;
    
    
    % Input-Output function
    f_i_mc_uncertainty = heaviside(i_uncertainty_total) *  i_uncertainty_total;
    f_i_mc_inhibitory = heaviside(i_inhibitory_total) *  i_inhibitory_total;
    
    % Dynamical equations.
    y_mc_uncertainty(t+1) = y_mc_uncertainty(t) + (dt/tau_mc_uncertainty) * ...
        (-y_mc_uncertainty(t) +f_i_mc_uncertainty);
    y_mc_inhibitory(t+1) = y_mc_inhibitory(t) + (dt/tau_mc_inhibitory) * ...
        (-y_mc_inhibitory(t) +f_i_mc_inhibitory);
    
    overall_mc_input(t+1) = i_uncertainty_total;
    %%%%%%%%%%%%%%%%%%%%%%%%%%% UNCERTAINTY MODULE INTEGRATE END
end
return;
end % Function end

function [y_5,y_6,response_time] =  hand_module_integrate(y_7,y_8,response_time)

global dt; global trial_length;

if(response_time==0)
    y_6=[];
    y_5=[];
    return;
end


j_dec_hand = 1.5;     % Coupling strength from decision module to hand module
hand_inhibition = 2;% Inhibition strength
j_self_hand = 0;    % Self excitation strength
tau_hand = 50;      % Coupling constant for hand
gain_hand = 1;      % Gain of input-output function

% Initialise state vectors with zeros (for better performance)
y_5 = zeros(1,trial_length);
y_6 = zeros(1,trial_length);

% Single trial loop start
for t=1:trial_length-1
    
    inhibitory_control = (t<response_time) * 5000;
    
    % Total input coming in to left hand neural population
    i5_total = gain_hand * (-hand_inhibition * y_6(t) + j_self_hand ...
        * y_5(t) + j_dec_hand* y_7(t)) - inhibitory_control;

    f_i5_total = heaviside(i5_total) * i5_total;
    
    % Dynamical equation for left hand neural population
    y_5(t+1) = y_5(t) + ((dt/tau_hand) * (-y_5(t) + f_i5_total));
    
    
    % Total input coming in to left hand neural population
    i6_total = gain_hand* (-hand_inhibition * y_5(t)+ j_self_hand ...
        * y_6(t) + j_dec_hand * y_8(t)) - inhibitory_control;
    
    f_i6_total = heaviside(i6_total) * i6_total;
    
    % Dynamical equation for left hand neural population
    y_6(t+1) = y_6(t) + ((dt/tau_hand) * (-y_6(t) + f_i6_total));
    
end % Single trial loop end

return;
end % Function end

function [is_motor_correct, motor_decision_made, movement_duration] =  check_if_motor_target_reached(y_5,y_6,motor_target_threshold,response_time)

motor_decision_made = false;
is_motor_correct = false;
movement_duration = 0;

left_reached = find(y_5>=motor_target_threshold,1);
right_reached = find(y_6>=motor_target_threshold,1);

if(isempty(left_reached))
    left_reached=0;
end

if(isempty(right_reached))
    right_reached=0;
end

if(right_reached > left_reached)
    is_motor_correct = true;
    motor_decision_made=true;
    movement_duration = movement_duration+ right_reached-response_time;
elseif(left_reached > right_reached)
    motor_decision_made=true;
    movement_duration = movement_duration + left_reached-response_time;
end

return;
end


function [is_motor_com, is_late_com] = check_com(y_5,y_6,motor_decision_made,motor_target_threshold,is_motor_correct)
global trial_length;

is_motor_com = false;
is_late_com = false;
if(~motor_decision_made)
    return;
end
x_traj = y_5-y_6;

smoothed_trajectory = filter(ones(1,50)/50,1,x_traj);
delta_of_trajectory = diff(sign(smoothed_trajectory));

left_reached = find(y_5>=motor_target_threshold,1);
right_reached = find(y_6>=motor_target_threshold,1);

if(isempty(left_reached))
    left_reached=0;
end

if(isempty(right_reached))
    right_reached=0;
end

points_of_change = find(delta_of_trajectory);

if(size(points_of_change,2)<2)
    return;
end

point_of_change = points_of_change(end);
point_before_point_of_change=points_of_change(end-1);

for i = point_of_change:trial_length-1
    if(y_5(i)>=motor_target_threshold || y_6(i)>=motor_target_threshold)
        is_motor_com = true;
        break;
    end
end

%check if correct Change-of-Mind
if(is_motor_com && is_motor_correct && point_of_change > left_reached && left_reached ~=0)
    is_late_com = true;
end

%check if incorrect Change-of-Mind
if(is_motor_com && ~is_motor_correct && point_of_change > right_reached && right_reached ~=0)
    is_late_com = true;
end

end


%%%%%%
%%%analysis:
function [p_com,p_correct,p_com_correct,p_com_incorrect,p_late_com_correct,p_late_com_incorrect] = calculate_accuracy_and_pcom(dynamics_and_results)



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



coupled = isfield(dynamics_and_results,'is_first_correct');

for i=1:n_of_coherences
    coherence= coherences(i);
    
    [noncom_correct_counter, noncom_incorrect_counter,...
        com_correct_counter, com_incorrect_counter, com_late_correct_counter, com_late_incorrect_counter, nondecision_counter] = ...
        calculate_accuracies(dynamics_and_results, coherence, coupled);
    
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
    
end

early_com_correct = com_correct - late_com_correct;
early_com_incorrect = com_incorrect - late_com_incorrect;

p_early_com_correct = early_com_correct./all_trials;
p_early_com_incorrect = early_com_incorrect./all_trials;


return;

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
    calculate_accuracies(dynamics_and_results, coherence, coupled)

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

function mean_rt = calculate_mean_rt(dynamics_and_results)

coherences = get_coherence_levels(dynamics_and_results);
n_of_coherences= size(coherences,2);



initiation_time_mean_gather = zeros(1,n_of_coherences);
initiation_time_std_gather = zeros(1,n_of_coherences);


for i=1:n_of_coherences
    coherence= coherences(i);
    

    [initiation_time_mean_gather(i),initiation_time_std_gather(i),~,~] = ...
        calculate_mean_times(dynamics_and_results, coherence); 
end
%%%%%%%%
mean_rt = initiation_time_mean_gather;
end

function[initiation_time_mean, initiation_time_std, ...
    response_time_mean, response_time_std] = ...
    calculate_mean_times(dynamics_and_results, coherence)

coherences = get_coherence_levels(dynamics_and_results);

if(~ismember(coherence,coherences))
    error('coherence level not found in data');
end

counter = 0;


for i = 1:size(dynamics_and_results,1)
    if(dynamics_and_results(i).coherence_level ==...
            coherence && dynamics_and_results(i).motor_decision_made)
        counter = counter+1;
    end
end


initiation_time_gather = zeros(1,counter);
response_time_gather = zeros(1,counter);

j=1;

for i = 1:size(dynamics_and_results,1)
    if(dynamics_and_results(i).coherence_level == ...
            coherence && dynamics_and_results(i).motor_decision_made)
        initiation_time_gather(j)=dynamics_and_results(i).response_time -dynamics_and_results(i).stim_onset ;
        response_time_gather(j)=dynamics_and_results(i).response_time + dynamics_and_results(i).movement_duration;
        j=j+1;
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
%%%%%%