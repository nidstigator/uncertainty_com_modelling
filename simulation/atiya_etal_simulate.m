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


dt = 0.5;
trial_length = 4000/dt;
stim_onset = 900;
stim_duration = 850;


%%%%Threshold values:
decision_module_threshold = 35.5;
motor_target_threshold = 17.4;
%%%%%%%%%%%%%%%%%%%%%%


prompt = 'Please enter the coherence level, or please enter 100 for all: \n';
coherences = input(prompt);

clearvars -except coherences motor_target_threshold decision_module_threshold ...
    stim_duration stim_onset trial_length dt 

close all;

coherence_string = [num2str(coherences) '_coh'];

if(coherences==100)
    coherences = [0,3.2,6.4,12.8,25.6,51.2];
    coherence_string = 'all_coh';
end

prompt = 'Please enter the number of trials: \n';
n_of_trials = input(prompt);

tStart=tic;


% This is the struct where all the results are saved.
% The struct is initiated before the beginning of any simulation.
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


across_trial_string = 'normal';
trial_type_string = 'reaction_time';

filename = [ 'data' '_' coherence_string '_' trial_type_string '_' across_trial_string '_'  num2str(n_of_trials) '_' strrep(strrep(datestr(now, 'dd-mm-yyyy HH-MM-SS'),' ', '_'), '-','_') '.mat' ];



index = 0;
no_of_blocks = size(coherences,2);

% Main nested for loop.
% Each evidence quality is assumed to be a block of trials.
for k = 1: no_of_blocks
    
    coh = coherences(k);
    
    % Trials for loop.
    for i = 1:n_of_trials
        index = index+1;
        
        %This function integrates the sensorimotor module and uncertainty
        %module.
        [y_1,y_2,s1,s2,y_mc_hu,y_mc_lu,decision_module_crossed,...
            response_time,overall_mc_input] = nonlinear_decision_module_and_mc_integrate...
            (coh,decision_module_threshold);
        
        %This function integrates the motor module.
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
            
            %This function checks for the result of the trial.
            [is_motor_correct, motor_decision_made, movement_duration] = ...
                check_if_motor_target_reached...
                (y_5,y_6,motor_target_threshold,response_time);
            
            %This function checks for change-of-mind.
            [is_motor_com,is_late_com] = check_com(y_5,y_6,motor_decision_made,motor_target_threshold,is_motor_correct);
        end
        
       
        %%%% Save trial in struct:
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
        clearvars -except dynamics_and_results k i dt trial_length stim_onset...
            stim_duration  decision_module_threshold ...
            motor_target_threshold coh n_of_trials tStart filename ...
            coherences index  no_of_blocks
        %%%%
    end
    
end
clearvars -except dynamics_and_results tStart filename

tElapsed=toc(tStart);

fprintf('Time taken to run this: %f \n',tElapsed)

clearvars tStart tElapsed;

function [y_1,y_2,s1,s2,y_mc_uncertainty,y_mc_inhibitory,decision_module_crossed,response_time,overall_mc_input] = nonlinear_decision_module_and_mc_integrate(coh,decision_module_threshold)

global dt; global trial_length; global stim_onset;

% Initialise response time.
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

gain_uncertainty = 1;   % Gain of input-output function (cf. Wong & Wang, 2006)
tau_mc_uncertainty = 150;        % Membrane time constant (cf. Wilson and Cowan, 1972)
gain_inh = 1;           % Gain of input-output function (cf. Wong & Wang, 2006)
tau_mc_inhibitory = 150;        % Membrane time constant (cf. Wilson and Cowan, 1972)

j_mc0 = 0.009;          % Excitation constant from uncertainty to sensorimotor
j_mc_v_inh = 1;         % Coupling strength from decision module to inhibtory neural population
mu= 30;                 % Tonic constant bias input used for uncertainty-encoding I-O function.
j_u_inh = 0.5;          % Inhibition strength constant from inhibitory neural population to to uncertainty-encoding neural population
j_self_uncertainty = 0;          % Self excitation constant for uncertainty populations.

t_u = 500;              % Integration delay after simtulus onset for uncertinaty-encoding population
t_inh = 500;            % Integration delay after simtulus onset for Inhibtory population
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

% Booleans to indicate decision threshold crossing in sensorimotor
% module.
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
    JN11 = 0.2440182353; JN22 = 0.2440182353; % 1.59
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
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%% UNCERTAINTY MODULE INTEGRATE:
    % Conditional inhibitory control.
    inh_top_down_inhibition = ((stim_onset+t_u)/dt>t) * 1000;
    uncertainty_top_down_inhibition = ((stim_onset+t_inh)/dt>t) * 1000;
    
    custom_top_down_inhibition = (decision_module_crossed) * 3000;
    
    % Total input coming in to uncertainty population (inhibitory and uncertainty)
    i_uncertainty_total = gain_uncertainty * (mu - ...
        j_u_inh * y_mc_inhibitory(t)) - inh_top_down_inhibition - custom_top_down_inhibition;
    
    i_inhibitory_total = gain_inh * (j_mc_v_inh * (y_1(t) +...
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
    %%%%%%%%%%%%%%%%%%%%%%%%%%%UNCERTAINTY MODULE INTEGRATE END
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
hand_inhibition = 2;  % Inhibition strength
j_self_hand = 0;      % Self excitation strength
tau_hand = 50;        % Coupling constant for hand
gain_hand = 1;        % Gain of input-output function

% Initialise state vectors with zeros (for better performance)
y_5 = zeros(1,trial_length);
y_6 = zeros(1,trial_length);

% Single trial loop start
for t=1:trial_length-1
    
    %TODO:Check this
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

% smoothed_trajectory = x_traj;
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


%take the last two points of change. 
point_of_change = points_of_change(end);

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