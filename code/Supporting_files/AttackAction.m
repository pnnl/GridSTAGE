% This code implements the attacks: Ramp, Step and Poisoning
% This code is called every time step during the course of the attack
% At each time step, the actual Voltage magnitude, voltage angle, frequency
% and rate of change of frequency measurements are corrupted according to
% the user-defined attack characteristics. 

if strcmp(AT,'Poisoning')
    
    if t > attack.start_time_in_sec && t < attack.end_time_in_sec
        row_location = find(PMU_samples == k, 1);
        row_location_fdot = find(PMU_samples_fdot == k, 1);
        
        if ~isempty(row_location)
            PMU.Vm(row_location, AttackLocation) = PMU.Vm(row_location,AttackLocation)+(AttackVector.MultiplicativeFactor(:,k-AttackVector.attack_start_step))';
            PMU.Va(row_location, AttackLocation) = PMU.Va(row_location,AttackLocation)+(AttackVector.MultiplicativeFactor(:,k-AttackVector.attack_start_step))';
            PMU.f(row_location, AttackLocation)  = PMU.f(row_location,AttackLocation)+(AttackVector.MultiplicativeFactor(:,k-AttackVector.attack_start_step))';
        end
        
        if ~isempty(row_location_fdot)
            PMU.fdot(row_location_fdot, AttackLocation)  = PMU.fdot(row_location_fdot,AttackLocation)+(AttackVector.MultiplicativeFactor(:,k-AttackVector.attack_start_step))';
        end
    end
    
elseif strcmp(AT, 'Ramp') || strcmp(AT, 'Step')
    % Ramp or Step attacks 
    if t > attack.start_time_in_sec && t < attack.end_time_in_sec
        row_location = find(PMU_samples == k, 1);
        row_location_fdot = find(PMU_samples_fdot == k, 1);
        
        if ~isempty(row_location)
            PMU.Vm(row_location, AttackLocation) = PMU.Vm(row_location,AttackLocation).*((AttackVector.MultiplicativeFactor(:,k-AttackVector.attack_start_step))');
            PMU.Va(row_location, AttackLocation) = PMU.Va(row_location,AttackLocation).*((AttackVector.MultiplicativeFactor(:,k-AttackVector.attack_start_step))');                                    
            PMU.f(row_location, AttackLocation)  = PMU.f(row_location,AttackLocation).*((AttackVector.MultiplicativeFactor(:,k-AttackVector.attack_start_step))');
        end
        
        if ~isempty(row_location_fdot)
            PMU.fdot(row_location_fdot, AttackLocation)  = PMU.fdot(row_location_fdot,AttackLocation).*AttackVector.MultiplicativeFactor(k-AttackVector.attack_start_step);
        end
    end
    
elseif strcmp(AT, 'Trapezoid')  
    % Trapezoidal attack
    
    % Increasing ramp    
    if t > attack.start_time_in_sec && t < attack.start_time_in_sec + Trapezoid.initial_slope 
        row_location = find(PMU_samples == k, 1);
        row_location_fdot = find(PMU_samples_fdot == k, 1);
        
        if ~isempty(row_location)
            PMU.Vm(row_location, AttackLocation) = PMU.Vm(row_location,AttackLocation).*(AttackVector.IncreasingRampMultiplicativeFactor(:,k-AttackVector.increasing_ramp_attack_start_step))';
            PMU.Va(row_location, AttackLocation) = PMU.Va(row_location,AttackLocation).*(0.2+(AttackVector.IncreasingRampMultiplicativeFactor(:,k-AttackVector.increasing_ramp_attack_start_step))');
            % PMU.f(row_location, AttackLocation)  = PMU.f(row_location,AttackLocation).*(AttackVector.IncreasingRampMultiplicativeFactor(:,k-AttackVector.increasing_ramp_attack_start_step))';
        end
        
        if ~isempty(row_location_fdot)
            PMU.fdot(row_location_fdot, AttackLocation)  = PMU.fdot(row_location_fdot,AttackLocation).*AttackVector.IncreasingRampMultiplicativeFactor(k-AttackVector.increasing_ramp_attack_start_step);
        end
    end
    
    % Step  
    if t >= attack.start_time_in_sec + Trapezoid.initial_slope  && t < attack.start_time_in_sec + Trapezoid.initial_slope + Trapezoid.intermediate_step
        row_location = find(PMU_samples == k, 1);
        row_location_fdot = find(PMU_samples_fdot == k, 1);
        
        if ~isempty(row_location)
            PMU.Vm(row_location, AttackLocation) = PMU.Vm(row_location,AttackLocation).*(AttackVector.IntermediateStepMultiplicativeFactor(:,k-AttackVector.intermediate_step_attack_start_step))';
            PMU.Va(row_location, AttackLocation) = PMU.Va(row_location,AttackLocation).*(0.2+(AttackVector.IntermediateStepMultiplicativeFactor(:,k-AttackVector.intermediate_step_attack_start_step))');
            % PMU.f(row_location, AttackLocation)  = PMU.f(row_location,AttackLocation).*(AttackVector.IntermediateStepMultiplicativeFactor(:,k-AttackVector.intermediate_step_attack_start_step))';
        end
        
        if ~isempty(row_location_fdot)
            PMU.fdot(row_location_fdot, AttackLocation)  = PMU.fdot(row_location_fdot,AttackLocation).*AttackVector.IntermediateStepMultiplicativeFactor(k-AttackVector.intermediate_step_attack_start_step);
        end
    end
    
    % Decreasing ramp    
    if t >= attack.start_time_in_sec + Trapezoid.initial_slope + Trapezoid.intermediate_step && t < attack.start_time_in_sec + Trapezoid.initial_slope + Trapezoid.intermediate_step + Trapezoid.final_slope
        row_location = find(PMU_samples == k, 1);
        row_location_fdot = find(PMU_samples_fdot == k, 1);
        
        if ~isempty(row_location)
            PMU.Vm(row_location, AttackLocation) = PMU.Vm(row_location,AttackLocation).*(AttackVector.DecreasingRampMultiplicativeFactor(:,k-AttackVector.decreasing_ramp_attack_start_step))';
            PMU.Va(row_location, AttackLocation) = PMU.Va(row_location,AttackLocation).*(0.2+(AttackVector.DecreasingRampMultiplicativeFactor(:,k-AttackVector.decreasing_ramp_attack_start_step))');
            % PMU.f(row_location, AttackLocation)  = PMU.f(row_location,AttackLocation).*(AttackVector.DecreasingRampMultiplicativeFactor(:,k-AttackVector.decreasing_ramp_attack_start_step))';
        end
        
        if ~isempty(row_location_fdot)
            PMU.fdot(row_location_fdot, AttackLocation)  = PMU.fdot(row_location_fdot,AttackLocation).*AttackVector.DecreasingRampMultiplicativeFactor(k-AttackVector.decreasing_ramp_attack_start_step);
        end
    end
    
elseif strcmp(AT, 'Freezing')
    % Freezing attack  
    if t > attack.start_time_in_sec && t < attack.end_time_in_sec
        row_location = find(PMU_samples == k, 1);
        row_location_fdot = find(PMU_samples_fdot == k, 1);
        
        if ~isempty(row_location)
            Freezing.time_point = 1800;
            
            check_frequency_freeze = find(Freezing.states == 'Frequency');
            check_voltage_freeze   = find(Freezing.states == 'Voltage');
            check_angle_freeze     = find(Freezing.states == 'Angle');
            
            if check_voltage_freeze
                PMU.Vm(row_location, AttackLocation) = PMU.Vm(Freezing.time_point,AttackLocation); 
            end
            if check_angle_freeze
                PMU.Va(row_location, AttackLocation) = PMU.Va(Freezing.time_point,AttackLocation); 
            end
            if check_frequency_freeze
                PMU.f(row_location, AttackLocation)  = PMU.f(Freezing.time_point,AttackLocation);
            end
                        
        end        
    end
    
elseif strcmp(AT, 'PacketDrop')
    row_location = find(PMU_samples == k, 1);
    PMU.Vm(row_location, AttackLocation) = PMU.Vm(row_location,AttackLocation).*BernoulliProcess(row_location);
    PMU.Va(row_location, AttackLocation) = PMU.Va(row_location,AttackLocation).*BernoulliProcess(row_location);
    row_location = find(PMU_samples_f == k, 1);
    PMU.f(row_location, AttackLocation)  = PMU.f(row_location,AttackLocation).*BernoulliProcess(row_location);
    row_location = find(PMU_samples_fdot == k, 1);
    PMU.fdot(row_location, AttackLocation)  = PMU.fdot(row_location,AttackLocation).*BernoulliProcess(row_location);
    
elseif strcmp(AT, 'Latency')
    row_location = find(PMU_samples == k, 1);
    PMU.Vm(NewTimes(row_location), AttackLocation) = PMU.Vm(row_location,AttackLocation);
    PMU.Vm(row_location:NewTimes(row_location),AttackLocation) = zeros(length(row_location:NewTimes(row_location)), length(AttackLocation));
            
end