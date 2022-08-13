% Author: Sai Pushpak Nandanoori (uses PST)
% Date created: October 22, 2019
% Updated: April 8, 2020 by Seemita Pal
% updated: April 27, 2020 by Sai Pushpak
% updated: May 12, 2020 by Sai Pushpak (added trapezoidal attacks
% and freezing attacks)
% updated: November 5, 2020 by Sai Pushpak (generating transient time-series
% data by making load changes and undoing load changes - basically the idea
% is to achieve the pre-disturbance equilibrium for post-disturbance)
%--------------------------------------------------------------------------
% Main code follows from here
clear all; clear global;
close all; % clc;
% PST_path = genpath('Matlab_PST'); % Provide filepath to Matlab PST folders
% addpath(PST_path)
SF_path = genpath('Supporting_files');
addpath(SF_path)
pst_path = genpath('/Users/nand311/OneDrive - PNNL/Documents/MATLAB/PST');
addpath(pst_path)
rng(2343,'v5uniform') % The seeding can be changed for different batch runs
% =========================================================================
% For the load modulation case (Laurentiu Dan Marinovici - 2017/02/14)
global h_sol TimeStep_of_simulation load_change_parameters n_lhc_samples n_lc_events_per_scenario
global fmeas_con busDeltaAngle busDeltaFreqPU dBusDeltaFreqPU busFreq dbusFreq
global line amount_of_load_change_per_scen loads_undergoing_change_sorted
% -------------------------------------------------------------------------
% For AGC
global agc_ratio agc_zm freq_zm agc_z1 agc_z2 agc_zm1 agc_zm2 agc_z
global agc_time_step agc_control gen_inputs load_changes num_area
% -------------------------------------------------------------------------
global ace1 ace2 S1_agcm S2_agcm area1_buses area2_buses
global basmva simParams
global PMU_locations PMU_samples PMU_samples_f PMU_samples_fdot
global AttackLocation AttackVector Freezing
global PMU SCADA n_PMU SCADA_row_location PMU_SCADA_Difference
global SCADA_locations SCADA_samples gen_changes
global PMU_area1_locations PMU_area2_locations n_loads_to_change
global ACE_data count count2 Trapezoid TieLineScheduledPowers
global PMU_SamplingFreq attack AttackTypes PMU_attack AT BernoulliProcess NewTimes
% =========================================================================

% User-defined parameters:
SavePlots = 1; % '1' for saving plots; '0' for not saving
% -------------------------------------------------------------------------
% System configurations
Network = '145'; % Choose the IEEE bus system: '9', '39', '68', '145'
% Load topology data to pick strategic load locations and generate transient
% data for power network
load(sprintf('TopologyData%s.mat',Network));
G = graph(A);

% -------------------------------------------------------------------------
identify_strategic_load_locations
% this identifies strategic load locations where the user can chose to
% introduce the disturbance
%{
buses_1_hop
buses_2_hop
buses_3_hop
buses_4_hop
load_buses_near_high_MW_gens
load_buses_near_low_MW_gens
load_buses_near_high_inertia_gens
load_buses_near_low_inertia_gens
%}
% -------------------------------------------------------------------------
agc_control = 0; % '1' enables AGC control; '0' disables AGC control
% AGC controls tuned only for IEEE 68 bus system
if ~strcmp(Network,'68'), agc_control = 0; end
agc_time_step = 0; % time interval in seconds between agc control
pss_control = 0; % '1' enables PSS control; '0' disables PSS control
num_area = 2; % 1- one area and 2 - two area
load_changes = 1; % '1' enables load changes; '0' disables load changes
gen_changes  = 0; % '1' enables gen changes;  '0' disables gen changes
TimeStep_of_simulation = 0.01; % in seconds
SimulationTime = 30; % in seconds
PMU_SamplingFreq  = 50; % Measurements every second
PMU_attack = 0;
% Number of fault scenarios (data with respect to each fault location will be saved as a single scenario)
n_fault_scenarios = 1;
% Number of fault-types (data with respect to each fault type will be saved as a single scenario)
n_fault_type = 1; % if there is no fault in the simulation, then chose n_fault_type = 1
FaultCases = {'three phase fault', 'line to ground fault',...
    'line-to-line to ground fault','line-to-line fault',...
    'loss of line with no fault', 'loss of load at bus','no fault'};

% NOTE: If a particular fault type needs to be simulated at a particular
% location, make changes accordingly to the sw_con matrix manually and
% chose, n_fault_scenarios = 1; n_fault_type = 1;
% -------------------------------------------------------------------------
% Number of load change scenarios
n_lc_scenarios = 1; % Number of load changes (== # num of scenarios corresponding to the load changes)
lhs_sd = [5 9 13];
%--------------------------------------------------------------------------
Initialization
%--------------------------------------------------------------------------
loads_undergoing_change = load_buses_near_high_MW_gens;
n_loads_to_change = length(loads_undergoing_change);
tmp_lvs = nominal_load_values(loads_undergoing_change);
loads_undergoing_change_sorted = zeros(n_loads_to_change,1);
for i_l = 1:n_loads_to_change
    loads_undergoing_change_sorted(i_l) = find(load_locations == loads_undergoing_change(i_l));
end

for i_lhs_sd = 1%:length(lhs_sd)
    % load changes -- based on Latin Hyper Cube sampling
    amount_of_load_changes = (lhsnorm(tmp_lvs', diag(abs(tmp_lvs)*lhs_sd(i_lhs_sd)),5000))'...
        - nominal_load_values(loads_undergoing_change);
    %{
    temp = amount_of_load_changes(1,:);
    [~, indx] = sort(abs(temp),'descend');
    amount_of_load_change = temp(:,indx(randi(300,n_lc_scenarios,1)));
    %}
    %
    amount_of_load_change = zeros(size(tmp_lvs,1),n_lc_scenarios);
    for i = 1:size(tmp_lvs,1)
        temp = amount_of_load_changes(i,:);
        [~, indx] = sort(abs(temp),'descend');
        amount_of_load_change(i,:) = temp(:,indx(randi(100,n_lc_scenarios,1)));
    end
    %}
    %--------------------------------------------------------------------------

    for i_load_changes = 1:n_lc_scenarios
        % The below load_change_parameters start_time and end_time properties can be either
        % scalar or vector [depending on how many load changes are needed in one run].
        % %% start_time is chosen randomly
        % load_change_parameters.start_time = diag(sort(randi(SimulationTime-15,n_lc_events),'ascend'));
        % %% start_time is chosen deterministically
        load_change_parameters.start_time = 1;
        load_change_parameters.end_time   = load_change_parameters.start_time + 0.25; % []

        amount_of_load_change_per_scen = amount_of_load_change(:,i_load_changes);

        for i_fault_location = 1:n_fault_scenarios
            for i_fault_type = 1:n_fault_type

                %{
                row 1 col1  simulation start time (s) (cols 2 to 6 zeros)
                      col7  initial time step (s)
                row 2 col1  fault application time (s)
                      col2  bus number at which fault is applied
                      col3  bus number defining far end of faulted line
                      col4  zero sequence impedance in pu on system base
                      col5  negative sequence impedance in pu on system base
                      col6  type of fault  - 0 three phase
                                           - 1 line to ground
                                           - 2 line-to-line to ground
                                           - 3 line-to-line
                                           - 4 loss of line with no fault
                                           - 5 loss of load at bus
                                           - 6 no action
                      col7  time step for fault period (s)
                row 3 col1  near end fault clearing time (s) (cols 2 to 6 zeros)
                      col7  time step for second part of fault (s)
                row 4 col1  far end fault clearing time (s) (cols 2 to 6 zeros)
                      col7  time step for fault cleared simulation (s)
                row 5 col1  time to change step length (s)
                      col7  time step (s)
                %}
                %{
                sw_con = [...
                    0.0    0    0    0    0    0      0.005; %sets intitial time step
                    0.05  line(i_fault_location,1)   line(i_fault_location,2)   0    0   i_fault_type       0.005; %apply three phase fault at bus 1, on line 1-2
                    0.1 0    0    0    0    0      0.005; %clear fault at bus
                    0.2 0    0    0    0    0      0.005; %clear remote end
                    0.5 0    0    0    0    0      0.005; % increase time step
                    %             1.0  0    0    0    0    0    0.005; % increase time step
                    SimulationTime  0    0    0    0    0    0]; % end simulation
                %}
                sw_con = [...
                    0.0   0    0    0    0    0      0.005; %sets intitial time step
                    0.05  6    7    0    0    6      0.005; %apply three phase fault at bus 1, on line 1-2
                    0.1   0    0    0    0    0      0.005; %clear fault at bus
                    0.2   0    0    0    0    0      0.005; %clear remote end
                    0.5   0    0    0    0    0      0.005; % increase time step
                    %             1.0  0    0    0    0    0    0.005; % increase time step
                    SimulationTime  0    0    0    0    0    0]; % end simulation
                % Changing the time-step of simulation to user-defined value
                sw_con(:,end) = TimeStep_of_simulation;

                FaultType = FaultCases(sw_con(2,6)+1);
                fprintf('[INFO] Fault: %s\n',upper(FaultType{1}))

                % Saving the Scenario description.txt file based on the user inputs:
                Bus1      = sw_con(2,2);
                Bus2      = sw_con(2,3);

                % Changing the time-step of simulation to user-defined value
                sw_con(:,end) = TimeStep_of_simulation;
                %--------------------------------------------------------------------------
                scenIdx_local = 0;

                tic
                if num_area == 1
                    ACE_data.k = []; ACE_data.myACE = [];
                elseif num_area == 2
                    ACE_data.tielinedifferences = [];
                    ACE_data.myACE1 = []; ACE_data.myACE2 = [];  ACE_data.k = [];
                end
                scenIdx_local = scenIdx_local+1;
                %----------------------------------------------------------
                scenario_description
                %------------------------------------------------------------------------
                %         load('scheduled_tie_powers.mat')
                % assuming [Buses 1-2; Buses 1-27; Buses 9-8]
                % for IEEE 68 bus system
                %         TieLineScheduledPowers  = [P_1_2_scheduled; P_1_27_scheduled; -P_8_9_scheduled];
                %----------------------------------------------------------
                s_simu_PD % Main code running power flow and time-domain simulations
                %----------------------------------------------------------
                data_processing_and_plotting
                %----------------------------------------------------------
                % Save the plots and data
                if SavePlots
                    FigureFile = sprintf('%s/PMU_V',scenDir);
                    saveas(voltage_plot,FigureFile,'jpg')
                    FigureFile = sprintf('%s/PMU_F',scenDir);
                    saveas(frequency_plot,FigureFile,'jpg')
                end
                close all;
                save(sprintf('%s/PMUData.mat',scenDir), 'lmod_sig','lmod_st','loads_undergoing_change_sorted','tmp_lvs','PMU','PMU_SamplingFreq','SimulationTime', 'fmeas_con', 'TimeStep_of_simulation', 'AttackLocation', 't', 'PMU_samples', 'pelect');
                % save(sprintf('%s/SCADAData.mat',scenDir), 'SCADA','Required_NumofMeasEvery2Seconds_SCADA','SimulationTime','Vmeas_con');
                % save(sprintf('%s/ControlInputData.mat',scenDir),'tg_sig')
                et = toc;
                ets = num2str(et);
                disp(['elapsed time = ' ets 's'])
                fprintf('------------------------*---------------------------\n')
            end
        end
    end
end