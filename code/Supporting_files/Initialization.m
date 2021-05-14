% Printing statements describing scenario
if load_changes, disp('[INFO] Load changes: YES'); else, disp('[INFO] Load changes: NO'); end
if pss_control, disp('[INFO] Power system stabilizer: YES'); else, disp('[INFO] Power system stabilizer: NO'); end
if agc_control, disp('[INFO] Automatic generation control: YES'); else, disp('[INFO] Automatic generation control: NO'); end

if PMU_attack 
    disp('[INFO] PMU attack: YES'); 
    fprintf('[INFO] Type of attack: %s\n',upper(AT)); 
    % fprintf('[INFO] Attack start time in seconds: %0.2f\n', attack.start_time_in_sec)
    % fprintf('[INFO] Attack end time in seconds: %0.2f\n', (attack.start_time_in_sec+attack.duration_in_sec))
else 
    disp('[INFO] PMU attack: NO'); 
end

% =========================================================================
switch Network
    case '9'
        dfile = 'd3m9bm.m';
%         dfile = 'd9bus.m';
    case '39'
        dfile = 'datane.m';
    case '68'
        dfile = 'data16m_pd.m';
        % dfile = 'data16m_original.m';
    case '145'
        dfile = 'data50m.m';
end
%--------------------------------------------------------------------------
pst_var % set up global variables 
jay = sqrt(-1);
%--------------------------------------------------------------------------
% Setting up for use in "s_simu_PD.m"
draw_voltage_plot = 0; % Plot voltage in real time
svc_dc=[];tcsc_dc=[];dcr_dc=[];dci_dc=[];
set(0,'DefaultAxesFontSize',14)
set(0,'DefaultTextInterpreter','Latex')
% =========================================================================
pathname = pwd;
if pathname == 0
   error('you must select a valid data file')
else
   lfile =length(dfile);
   % strip off .m and convert to lower case
   dfile = lower(dfile(1:lfile-2));
   eval(dfile);   
end

% PMU location information:
PMU_locations = fmeas_con(:,2)';
n_PMU = length(PMU_locations);

% SCADA location information:
SCADA_locations = Vmeas_con(:,2);

load_locations = find(bus(:,6)>0);
num_loads = length(load_locations);

nominal_load_values = bus(:,6); % real-power [load values] 

% Generator input signals
gen_inputs = [];

n_lhc_samples = n_lc_scenarios*n_lc_events_per_scenario;    
mac_con(:,19) = 3; % Choose the model of generator (Synchronous machine)
%------------------------------------------------------------------------
% Defining sampling rates for PMU and SCADA measurements:
% PMU 
NumberofMeasPerSecond               = 1/TimeStep_of_simulation;
Number_of_TimeSteps                 = SimulationTime/TimeStep_of_simulation;
Required_NumberofMeasPerSecond_PMU_given = PMU_SamplingFreq;

% Target PMU sampling rate based on current time step of the simulation. 
Tar_PMU_sr                          = (1/PMU_SamplingFreq)/TimeStep_of_simulation; 
myflag = 0;

while mod(Tar_PMU_sr,1)~=0
    myflag = 1;
    PMU_SamplingFreq = PMU_SamplingFreq+1;
    Tar_PMU_sr                          = (1/PMU_SamplingFreq)/TimeStep_of_simulation; 
end
if myflag
    disp('[INFO] Chosen number of PMU samples is not feasible. So, choosing a different number of samples per second.')
    fprintf('[INFO] Given number of PMU samples every second: %d\n', Required_NumberofMeasPerSecond_PMU_given);
    fprintf('[INFO] Modified number of PMU samples every second: %d\n', PMU_SamplingFreq);
    fprintf('------------------------*---------------------------\n')
end

%------------------------------------------------------------------------
% SCADA 
Required_NumofMeasEvery2Seconds_SCADA = 1; 
% Target SCADA sampling rate based on current time step of the simulation. 
Tar_SCADA_sr = (2/Required_NumofMeasEvery2Seconds_SCADA)/TimeStep_of_simulation;
while mod(Tar_SCADA_sr,1) ~= 0
    Required_NumofMeasEvery2Seconds_SCADA = Required_NumofMeasEvery2Seconds_SCADA+1;
    Tar_SCADA_sr = (2/Required_NumofMeasEvery2Seconds_SCADA)/TimeStep_of_simulation;
end
SCADA_sample_index = 0;

%------------------------------------------------------------------------
% Entries to be stored inside the CSV file per scenario
scenDes = cell(2,15);
scenDes{1,1}  = 'Fault';
scenDes{1,2}  = 'Attack';
scenDes{1,3}  = 'Load Changes';

scenDes{1,4}  = 'AGC';
scenDes{1,5}  = 'PSS';

scenDes{1,6}  = 'Fault location';
scenDes{1,7}  = 'Fault type';

scenDes{1,8}  = 'Attack location';
scenDes{1,9}  = 'Attack type';
scenDes{1,10} = 'Attack start time (sec)';
scenDes{1,11} = 'Attack end time (sec)';
scenDes{1,12} = 'Attack duration (sec)';
if strcmp(AT, 'Freezing')   
    scenDes{1,13} = 'Time point frozen'; 
else
    scenDes{1,13} = 'Attack amplitude (%)';
end

scenDes{1,14} = 'Location(s) of load changes';
scenDes{1,15} = 'Magnitude(s) of load changes';
scenDes{1,16} = 'Start time(s) for load changes';
scenDes{1,17} = 'End time(s) for load changes';

%------------------------------------------------------------------------
% Entries to be stored inside the CSV file for one batch run 
n_scenarios = n_attack_scenarios*n_attacks_on_magnitude*n_lc_scenarios*n_attacks_on_duration_of_attack*n_lc_scen;
scenDes_full = cell(n_scenarios+1,22);
scenDes_full{1,1} = 'Scenario Number';

scenDes_full{1,2} = scenDes{1,1};
scenDes_full{1,3} = scenDes{1,2};
scenDes_full{1,4} = scenDes{1,3};

scenDes_full{1,5} = scenDes{1,4};
scenDes_full{1,6} = scenDes{1,5};

scenDes_full{1,7} = scenDes{1,6};
scenDes_full{1,8} = scenDes{1,7};

scenDes_full{1,9} = scenDes{1,8};
scenDes_full{1,10} = scenDes{1,9};
scenDes_full{1,11} = scenDes{1,10};
scenDes_full{1,12} = scenDes{1,11};
scenDes_full{1,13} = scenDes{1,12};
scenDes_full{1,14} = scenDes{1,13};

scenDes_full{1,15} = scenDes{1,14};
scenDes_full{1,16} = scenDes{1,15};

scenDes_full{1,17} = scenDes{1,16};
scenDes_full{1,18} = scenDes{1,17};

%------------------------------------------------------------------------
% Compute P and Q from PMU data and save them into PMU structure
save_PMU_P_and_Q = 0;
% In this work, it is assumed that only PMUs are under attack and PMU, 
% SCADA measurements are available at every bus. This assumption is not 
% valid in general. Under this assumption, one can compute certain metrics 
% to identify the attacked PMUs. 
% Consider the difference between the PMU and SCADA datasets whereever 
% the SCADA measurements are available and take their various euclidean 
% norms such as 1-norm, 2-norm, Inf-norm. 
% voltage magnitude of each bus from the PMU data are compared against 
% SCADA voltage magnitude measurements. 
compute_metrics = 0;
if compute_metrics
    scenDes_full{1,19} = 'Most affected buses (1-norm)';
    scenDes_full{1,20} = 'Most affected buses (2-norm)';
    scenDes_full{1,21} = 'Most affected buses (Infty-norm)';
    scenDes_full{1,22} = 'Max impact locations (1-norm)';
    scenDes_full{1,23} = 'Max impact locations (2-norm)';
    scenDes_full{1,24} = 'Max impact locations (Infty-norm)';
end 