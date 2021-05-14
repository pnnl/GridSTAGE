% In this code, we define a multiplicative vector that is multiplied or
% added to the actual PMU data to mimic the attacks such as ramp/step or
% random poisoning

% Determining simulation steps when attack starts or ends
if strcmp(AT, 'Ramp') || strcmp(AT, 'Step') || strcmp(AT, 'Posioning') || strcmp(AT, 'Freezing')
    AttackVector.attack_start_step  = attack.start_time_in_sec/TimeStep_of_simulation;
    AttackVector.attack_end_step    = attack.end_time_in_sec/TimeStep_of_simulation;
elseif strcmp(AT,'Trapezoid')
    AttackVector.increasing_ramp_attack_start_step   = attack.start_time_in_sec/TimeStep_of_simulation;
    AttackVector.increasing_ramp_attack_end_step     = (attack.start_time_in_sec + Trapezoid.initial_slope)/TimeStep_of_simulation;
    AttackVector.intermediate_step_attack_start_step = AttackVector.increasing_ramp_attack_end_step + 1;
    AttackVector.intermediate_step_attack_end_step   = (attack.start_time_in_sec + Trapezoid.initial_slope + Trapezoid.intermediate_step)/TimeStep_of_simulation;
    AttackVector.decreasing_ramp_attack_start_step   = AttackVector.intermediate_step_attack_end_step + 1;
    AttackVector.decreasing_ramp_attack_end_step     = (attack.start_time_in_sec + Trapezoid.initial_slope + Trapezoid.intermediate_step + Trapezoid.final_slope)/TimeStep_of_simulation;
end
%------------------------------------------------------------------------
% number of attack locations
l_al = length(AttackLocation);
switch AT
    case 'Ramp'
        AttackVector.MultiplicativeFactor = [];
        % For the duration of the attack, creating a multiplicativeFactor
        % whose size is equal is to the number of attack locations x attack
        % duration (in time steps)
        % As this is a ramp attack, the slope of the attack is increased
        % during the course of the attack.
        for i_al = 1:l_al
%             AttackVector.MultiplicativeFactor = [AttackVector.MultiplicativeFactor; ...
%                 (linspace(1,attack.max_mod_frac(i_al),AttackVector.attack_end_step-AttackVector.attack_start_step))];
            AttackVector.MultiplicativeFactor = [AttackVector.MultiplicativeFactor; ...
                (linspace(1,attack.max_mod_frac(i_al),AttackVector.attack_end_step-AttackVector.attack_start_step))];
        end
    case 'Step'
        AttackVector.MultiplicativeFactor = [];
        % For the duration of the attack, creating a multiplicativeFactor
        % whose size is equal is to the number of attack locations x attack
        % duration (in time steps)
        % As this is a step attack, the magnitude of the attack is fixed
        % during the course of the attack.
        for i_al = 1:l_al
            AttackVector.MultiplicativeFactor = [AttackVector.MultiplicativeFactor; ...
                repmat(attack.max_mod_frac(i_al),1,AttackVector.attack_end_step-AttackVector.attack_start_step)];
        end
    case 'PacketDrop'
        % Stochastic process to create Packet drops
        Lambda = 0.98; % Probability of retaining the PMU data
        BernoulliProcess = binornd(1, Lambda*ones(SimulationTime*PMU_SamplingFreq,1));
    case 'Latency'
        % Stochastic process to mimic random time delays
        % deterministic time delay
        Lambda = 0.1; % Rate for Poisson distribution
        if strcmp(LP.type, 'Deterministic')
            NewTimes = (1:SimulationTime*PMU_SamplingFreq)'+LP.deterministic_time_shift;
        elseif strcmp(LP.type, 'Stochastic')
            DelayTimes = cumsum(poissrnd(Lambda,[SimulationTime*PMU_SamplingFreq,1]));
            NewTimes = (1:SimulationTime*PMU_SamplingFreq)' + DelayTimes;
        end
    case 'Poisoning'
        % For the duration of the attack, creating a multiplicativeFactor
        % whose size is equal is to the number of attack locations x attack
        % duration (in time steps)
        % As this is a poisoning attack, the multiplicativeFactor is added
        % to the actual data
        % As we don't want to make very abrupt changes due to poisoning
        % attack every second, the changes happen every AGC cycle.
        rep_frequency = 250;% agc_time_step*PMU_SamplingFreq;
        % Creating a random vector
        temp_randomf = randn(l_al,round(AttackVector.attack_end_step-AttackVector.attack_start_step)/rep_frequency);
        % Initializing a temporary attack vector
        temp_res_vec = zeros(l_al,AttackVector.attack_end_step-AttackVector.attack_start_step);
        for i_sub_attack_location = 1:l_al % corresponding to every attack location
            temp_random = temp_randomf(i_sub_attack_location,:);
            % Freezing the attack for the AGC Cycle
            temp_repeat  = repmat(temp_random',1,rep_frequency);
            % reshaping it to match with the duration of the attack
            temp_res_vec(i_sub_attack_location,:) = reshape(temp_repeat',1,AttackVector.attack_end_step-AttackVector.attack_start_step);
        end
        AttackVector.MultiplicativeFactor = poisoning.mean + poisoning.var.*temp_res_vec;
        %         AttackVector.MultiplicativeFactor = poisoning.mean + poisoning.var.*randn(length(AttackLocation),AttackVector.attack_end_step-AttackVector.attack_start_step);
    case 'Trapezoid'
        AttackVector.IncreasingRampMultiplicativeFactor    = [];
        AttackVector.IntermediateStepMultiplicativeFactor  = [];
        AttackVector.DecreasingRampMultiplicativeFactor    = [];
        % For the duration of the attack, creating a multiplicativeFactor
        % whose size is equal is to the number of attack locations x attack
        % duration (in time steps)
        
        % For increasing ramp attack, the slope of the attack is increased
        % during the course of the attack.
        for i_al = 1:l_al
            AttackVector.IncreasingRampMultiplicativeFactor = [AttackVector.IncreasingRampMultiplicativeFactor; ...
                (linspace(1,attack.initial_max_mod_frac(i_al),AttackVector.increasing_ramp_attack_end_step-AttackVector.increasing_ramp_attack_start_step))];
        end
        
        % For intermediate step attack, the slope of the attack is fixed
        % during the course of the attack.
        for i_al = 1:l_al
            AttackVector.IntermediateStepMultiplicativeFactor = [AttackVector.IntermediateStepMultiplicativeFactor; ...
                repmat(attack.intermediate_max_mod_frac(i_al),1,AttackVector.intermediate_step_attack_end_step-AttackVector.intermediate_step_attack_start_step)];
        end
        
        % For decreasing ramp attack, the slope of the attack is decreased
        % during the course of the attack.
        for i_al = 1:l_al
            AttackVector.DecreasingRampMultiplicativeFactor = [AttackVector.DecreasingRampMultiplicativeFactor; ...
                (linspace(attack.final_max_mod_frac(i_al),1,AttackVector.decreasing_ramp_attack_end_step-AttackVector.decreasing_ramp_attack_start_step))];
        end
        
end