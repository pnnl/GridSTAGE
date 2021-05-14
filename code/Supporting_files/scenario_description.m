Gen_Types = {'Electromechanical','Subtransient salient pole','Subtransient round rotor'}; 

% Creating Results and Scenario directory:
if ~exist('Results','dir'), mkdir('Results'); end
if PMU_attack 
    ResDir = sprintf('Results/IEEE%sbusSystem/%s',Network,AT);
else
    ResDir = sprintf('Results/IEEE%sbusSystem',Network);
end
if load_changes && ~PMU_attack
    ResDir = sprintf('Results/IEEE%sbusSystem/LoadChanges',Network);
end 
if ~exist(ResDir,'dir'), mkdir(ResDir); end

scenIdx = 1;
while exist(sprintf('%s/scenario%d',ResDir,scenIdx),'dir')
    scenIdx = scenIdx + 1;
end
scenDir = sprintf('%s/scenario%d',ResDir,scenIdx);
mkdir(scenDir)
% Saving the scenario information into the txt file
fileID = fopen(sprintf('%s/ScenarioDescription.txt',scenDir),'w'); 

fprintf(fileID, sprintf('IEEE %s bus system\n',Network));

fprintf(fileID, sprintf('Fault type: %s \n',FaultType{1}));
if ~strcmp(FaultType,'no fault')
    fprintf(fileID, sprintf('Fault location: between buses %d-%d \n \n', Bus1, Bus2));
end 
if load_changes 
    fprintf(fileID, 'Load changes: YES \n'); 
    fprintf(fileID, 'Start time of load changes: %0.2f \n', load_change_parameters.start_time);
    fprintf(fileID, 'End time of load changes: %0.2f \n', load_change_parameters.end_time);
    fprintf(fileID, 'Load buses where the load changes: %d \n', loads_undergoing_change);
    fprintf(fileID, 'Magnitude of load changes: %0.2f \n', amount_of_load_change_per_scen);    
else 
    fprintf(fileID, 'Load changes: NO \n'); 
end

if pss_control, fprintf(fileID, 'Power system stabilizer (PSS): YES \n'); else, fprintf(fileID, 'Power system stabilizer (PSS): NO \n'); end

if agc_control, fprintf(fileID, 'Automatic generation control (AGC): YES \n'); else, fprintf(fileID, 'Automatic generation control (AGC): NO \n'); end

if PMU_attack 
    fprintf(fileID, sprintf('PMU attack: YES\n')); 
    fprintf(fileID, sprintf('Type of attack: %s \n', AT)); 
    fprintf(fileID, sprintf('Attack Location: Buses %d, %d \n',AttackLocation)); 
    if ~strcmp(AT, 'Trapezoid')
        fprintf(fileID, sprintf('Start and end time of attack: %0.2f -- %0.2f seconds\n',attack.start_time_in_sec,attack.end_time_in_sec));   
    end 
    switch AT
        case 'Ramp'
            fprintf(fileID, sprintf('Attack magnitude (Percentage of deviation): %0.2f%% \n \n \n',(attack.max_mod_frac-1)*100));
        case 'Step' 
            fprintf(fileID, sprintf('Attack magnitude (Percentage of deviation): %0.2f%% \n \n \n',(attack.max_mod_frac-1)*100));
        case 'Posioning'
            fprintf(fileID, sprintf('Attack magnitude (variance of disturbance): %0.2f%% \n \n \n',data_poison.var));
            fprintf(fileID, sprintf('Attack magnitude (Mean of disturbance): %0.2f%% \n \n \n',data_poison.mean));
        case 'Trapezoid'
            fprintf(fileID, sprintf('Increasing ramp attack start and end time: %0.2f -- %0.2f seconds\n', attack.start_time_in_sec, attack.start_time_in_sec + Trapezoid.initial_slope));
            fprintf(fileID, sprintf('Step attack start and end time: %0.2f -- %0.2f seconds\n', attack.start_time_in_sec + Trapezoid.initial_slope, attack.start_time_in_sec + Trapezoid.initial_slope + Trapezoid.intermediate_step));
            fprintf(fileID, sprintf('Decreasing ramp attack start and end time: %0.2f -- %0.2f seconds\n', attack.start_time_in_sec + Trapezoid.initial_slope + Trapezoid.intermediate_step , attack.start_time_in_sec + Trapezoid.initial_slope + Trapezoid.intermediate_step + Trapezoid.final_slope));
    
            fprintf(fileID, sprintf('Increasing ramp attack magnitudes: %0.2f \n\n\n', Trapezoid.initial_attack_magnitude));
            fprintf(fileID, sprintf('Step attack magnitudes: %0.2f \n\n\n', Trapezoid.intermediate_attack_magnitude));
            fprintf(fileID, sprintf('Decreasing ramp attack magnitudes: %0.2f \n\n\n', Trapezoid.final_attack_magnitude));
        case 'Freezing'
            fprintf(fileID, sprintf('Frozen attack uses the state value at time point: %0.2f \n',(Freezing.time_point)*2*TimeStep_of_simulation));
            fprintf(fileID, sprintf('States frozen during the attack: %s \n', Freezing.states));
    end    
    
    fprintf(fileID, sprintf('\n \n \n'));
else
    fprintf(fileID, sprintf('PMU attack: NO\n')); 
end

fprintf(fileID, sprintf('Generator Type (at all generators): %s \n \n',Gen_Types{mac_con(1,19)}));

fprintf(fileID, 'Reproducing the sw_con variable information to understand the step time(s), simulation time(s), and event description. \n\n');
fprintf(fileID, 'sw_con = \n');
fprintf(fileID, '%0.2f %d %d %d %d %d %0.3f \n', sw_con');
fprintf(fileID, '\n \nPMU locations = \n');
fprintf(fileID, '%d \t', PMU_locations);
fprintf(fileID, '\n \nSCADA locations = \n');
fprintf(fileID, '%d \t', SCADA_locations);
fprintf(fileID, '\n \nNumber of PMU measurements per second) = \n');
fprintf(fileID, '%d \t', PMU_SamplingFreq);
fprintf(fileID, '\n \nNumber of SCADA measurements every 2 seconds) = \n');
fprintf(fileID, '%d \t', Required_NumofMeasEvery2Seconds_SCADA);
% fprintf(fileID, '\n \n \nNumber of time windows to save the simulated data:%d \n',Number_of_time_windows);
fclose(fileID);
% Saving the attack scenario description (ScenDes) as a csv file:
%{
scenDes = cell(2,13);
scenDes{1,1} = 'Normal';
scenDes{1,2} = 'Fault';
scenDes{1,3} = 'Attack presence';
scenDes{1,4} = 'Attack type';
scenDes{1,5} = 'Attack location';
scenDes{1,6} = 'AGC';
scenDes{1,7} = 'PSS';
scenDes{1,8} = 'Attack start time (sec)';
scenDes{1,9} = 'Attack duration (sec)';
scenDes{1,10} = 'Attack end time (sec)';
scenDes{1,11} = 'Attack amplitude (%)';
scenDes{1,12} = 'Fault location';
scenDes{1,13} = 'Fault type';
%}
scenDes{2,1} = ~strcmp(FaultType,'no fault');
scenDes{2,2} = PMU_attack;
scenDes{2,3} = load_changes;

scenDes{2,4} = agc_control;
scenDes{2,5} = pss_control;

if strcmp(FaultType,'no fault')
    scenDes{2,6} = 'None';
    scenDes{2,7} = 'None';
else
    scenDes{2,6} = string(sw_con(2,2:3));
    scenDes{2,7}  = FaultType{1};
end 

if PMU_attack     
    scenDes{2,8}  = num2str(AttackLocation); 
    scenDes{2,9}  = AT;
    scenDes{2,10} = attack.start_time_in_sec;
    scenDes{2,11} = attack.end_time_in_sec;
    scenDes{2,12} = attack.duration_in_sec;
    switch AT
        case 'Ramp'
            scenDes{2,13} = (attack.max_mod_frac-1)*100;    
        case 'Step' 
            scenDes{2,13} = (attack.max_mod_frac-1)*100;    
        case 'Posioning'
            scenDes{2,13} = data_poison.var;   
        case 'Freezing'
            scenDes{2,13} = (Freezing.time_point)*2*TimeStep_of_simulation; 
    end
else 
    scenDes{2,8}  = 'None';
    scenDes{2,9}  = 'None';
    scenDes{2,10} = 'None';
    scenDes{2,11} = 'None';
    scenDes{2,12} = 'None';
    scenDes{2,13} = 'None';
end

if load_changes 
    scenDes{2,14} = num2str(loads_undergoing_change');
    scenDes{2,15} = num2str(amount_of_load_change_per_scen');
    scenDes{2,16} = num2str(load_change_parameters.start_time);
    scenDes{2,17} = num2str(load_change_parameters.end_time);
else
    scenDes{2,14} = 'None';
    scenDes{2,15} = 'None';
    scenDes{2,16} = 'None';
    scenDes{2,17} = 'None';
end

scenDes_table = cell2table(scenDes);
 
% Write the table to a CSV file
writetable(scenDes_table,sprintf('%s/ScenarioDescription.csv',scenDir),'WriteVariableNames',0)
% writecell(scenDes,sprintf('%s/ScenarioDescription.csv',scenDir))
disp('[INFO] Scenario description saved.')
fprintf('Scenario: %d \n',scenIdx)                