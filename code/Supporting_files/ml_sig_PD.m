function f = ml_sig_PD(t,k)
% Syntax: f = ml_sig(t,k)
%4:40 PM 15/08/97
% defines modulation signal for lmod control
global lmod_sig n_lmod
% fprintf('(t,k) = (%0.3f,%d)\n',t,k);
% =====================================================================================
% For the load modulation case (Laurentiu Dan Marinovici - 2017/02/14)
% global line basmva load_con
global bus_int bus_v h_sol basrad load_change_parameters n_lc_events_per_scenario
global fmeas_con busDeltaAngle busDeltaFreqPU dBusDeltaFreqPU busFreq
global load_changes simParams PMU_SamplingFreq dbusFreq
global PMU_locations PMU_samples PMU_samples_f PMU_samples_fdot AttackTypes PMU_attack
global AttackLocation AttackVector PMU SCADA PMU_SCADA_Difference
global SCADA_locations SCADA_samples SCADA_row_location Freezing 
global attack AT BernoulliProcess NewTimes count2 Trapezoid n_loads_to_change
global amount_of_load_change_per_scen loads_undergoing_change_sorted 
% =====================================================================================

f = 0; % dummy variable
% fprintf('(t,k) = (%0.3f,%d)\n',t,k);
% =====================================================================================
% For the load modulation case (Laurentiu Dan Marinovici - 2017/02/14)
nominalFreq = basrad/2/pi;
measBusInd = bus_int(fmeas_con(:, 2)); % index of the HVDC line end buses
% busAngleVect = unwrap(angle(bus_v(measBusInd, 1:k)), [], 2);
busAngleVect = angle(bus_v(measBusInd, 1:k));
currentAngle = busAngleVect(:, end); % angle(bus_v(modBusInd, k));

if k == 1
    previousAngle = currentAngle;
else
    previousAngle = busAngleVect(:, end - 1); % angle(bus_v(modBusInd, k - 1));
end

busDeltaAngle(:, k) = currentAngle - previousAngle;
dBusDeltaFreqPU(:, k) = (-busDeltaFreqPU(:, k) + busDeltaAngle(:, k) / (h_sol * basrad)) ./ fmeas_con(:, 3);

% frequency at PMU buses, in the order of the buses listed in
% fmeas_con; this should be a 2 by k matrix, and the frequency difference
% at time k would be the difference between the row
busFreq(:, k) = (1 + busDeltaFreqPU(:, k)) .* (nominalFreq); % [bus 49; bus 24]

if k > 1
    dbusFreq(:, k-1) = ( busFreq(:,k) - busFreq(:,k-1) ) ./  (1/PMU_SamplingFreq) ;
end

% Finding the time instance at which the PMU data has to be saved
save_PMU_data = find(PMU_samples == k, 1);

% Saving the voltage magnitude and voltage angle data
if ~isempty(save_PMU_data)
    busVoltMag = abs(bus_v(PMU_locations,k)); % [pu]
    busVoltAng = angle(bus_v(PMU_locations,k))*180/pi; % degrees
    
    row_location = save_PMU_data;
    
    PMU.Vm(row_location,PMU_locations)   = busVoltMag';
    PMU.Va(row_location,PMU_locations)   = busVoltAng';
end

% finding the time instance at which the frequency has to be saved. Note
% that the voltage data and frequency data store at the same time instance
save_PMU_f_data = find(PMU_samples_f == k, 1);

% Saving the frequency data
if ~isempty(save_PMU_f_data)
    PMU.f(save_PMU_f_data, PMU_locations) = (busFreq(:,k))';
end

% Saving the rate of change of frequency
if k > 1
    save_PMU_fdot_data = find(PMU_samples_fdot == k, 1);
    if ~isempty(save_PMU_fdot_data)
        PMU.fdot(save_PMU_fdot_data, PMU_locations) = (dbusFreq(:,k-1))';
    end
end

% =========================================================================
% Finding the time instance at which the SCADA data has to be saved
save_SCADA_data = find(SCADA_samples == k, 1);

% Voltage magnitude and angle at SCADA buses at SCADA measuring instances
if ~isempty(save_SCADA_data)
    busVoltMag_S = abs(bus_v(SCADA_locations,k)); % [pu]
    SCADA_row_location = save_SCADA_data;
    SCADA.Vm(SCADA_row_location, SCADA_locations) = busVoltMag_S';
end
% =========================================================================
% Attack modeling for PMU data to explore effects on AGC control
if PMU_attack
    % disp('Attacking!')    
    AttackAction                    
end
% =========================================================================
% disp('Inside ml_sig!')
lmod_sig(:,k)=zeros(n_lmod,1);
% start_time = 15.0;
if n_lmod~=0
    lmod_sig(:,k)=zeros(n_lmod,1);
    %    lmod_sig(:,k) = 0.1*randn(n_lmod,1);
    if load_changes
%         disp('Inside Load Changes!')
%         for i_l = 1:n_loads_to_change  
%             fprintf('(t,k) = (%0.3f,%d)\n',t,k);
%             temp_indx = find(t >= load_change_parameters.start_time);
%             lmod_sig(loads_undergoing_change_sorted(temp_indx),k) = amount_of_load_change_per_scen(temp_indx);
            if t >= load_change_parameters.start_time
                lmod_sig(loads_undergoing_change_sorted,k) = amount_of_load_change_per_scen;
            end                         
            if ~isempty(load_change_parameters.end_time)
                if t <= load_change_parameters.end_time
                    lmod_sig(:,k)=zeros(n_lmod,1);
                end
            end       
%         end  
    end
end
% =========================================================================
% Time instant at which PMU measurement is computed
% Freeze the SCADA data from the previous time measurement

if ~isempty(save_PMU_data)
    PMU_SCADA_Difference(row_location,:) = PMU.Vm(row_location, PMU_locations) - SCADA.Vm(SCADA_row_location, PMU_locations);
end