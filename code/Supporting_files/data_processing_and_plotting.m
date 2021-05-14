% Processing data for PMUs and SCADA measurements and saving them to
% appropriate variables and plotting. 

%------------------------------------------------------------------------
% Processing data (PMU)
% For each node, identify the branches and find current through each branch
% to store them into PMUs
AllBranches = line(:,1:2);
for i_bus = 1:n_bus
    
    if sum(PMU_locations == i_bus)
        % Finding the number of branches connecting every node and
        % declaring the variables with appropriate sizes
        nodeBranches = [];
        nodeBranchIndx = [];
        for i_b = 1:length(AllBranches)
            if AllBranches(i_b,1) == i_bus || AllBranches(i_b,2) == i_bus
                nodeBranches = [nodeBranches; AllBranches(i_b,:)];
                nodeBranchIndx = [nodeBranchIndx; i_b];
            end
        end
        
        Im_node = zeros(length(t), length(nodeBranchIndx));
        Ia_node = zeros(length(t), length(nodeBranchIndx));
        I_description = cell(1, length(nodeBranchIndx));
        
        for i_b = 1:size(nodeBranches,1)
            if nodeBranches(i_b,1) ~= i_bus
                nodeBranches(i_b,:) = fliplr(nodeBranches(i_b,:));
                V2  = bus_v(bus_int(line(nodeBranchIndx(i_b),1)),:);
                V1  = bus_v(bus_int(line(nodeBranchIndx(i_b),2)),:);
            else
                V1  = bus_v(bus_int(line(nodeBranchIndx(i_b),1)),:);
                V2  = bus_v(bus_int(line(nodeBranchIndx(i_b),2)),:);
            end
            
            Sm_node = fromBusPowInj(nodeBranchIndx(i_b),:);
            I_node  = conj(Sm_node./V1);
            Im_node(:, i_b) = (abs(I_node))';
            Ia_node(:, i_b) = (angle(I_node))';
            I_description{i_b} = sprintf('I_%d_%d',nodeBranches(i_b,1),nodeBranches(i_b,2));
        end
        
        PMU.Im{i_bus} = Im_node(PMU_samples,:);
        PMU.Ia{i_bus} = Ia_node(PMU_samples,:);
        PMU.Id{i_bus} = I_description;
        
        % Additionally if the P and Q needs to be extracted from voltage
        % magnitude, voltage angle, current angle and current magnitudes. 
        if save_PMU_P_and_Q 
            jay = sqrt(-1);
            % ComplexCurrent = PMU.Im{i_bus} + jay*PMU.Ia{i_bus};
            [I_real, I_imag] = pol2cart(PMU.Ia{i_bus}, PMU.Im{i_bus});
            [V_real, V_imag] = pol2cart(deg2rad(PMU.Va(:,i_bus)),PMU.Vm(:,i_bus));
            ComplexCurrent = I_real + jay*I_imag;
            ComplexVoltage = V_real + jay*V_imag;
            ComplexPower   = ComplexVoltage.*conj(ComplexCurrent);

            PMU.P{i_bus} = real(ComplexPower);
            PMU.Q{i_bus} = imag(ComplexPower);
        end 
    else
        % This node doesn't have a PMU, so, no data
        PMU.Im{i_bus} = [];
        PMU.Ia{i_bus} = [];
        PMU.Id{i_bus} = [];
        if save_PMU_P_and_Q 
            PMU.P{i_bus}  = [];
            PMU.Q{i_bus}  = [];
        end 
    end
end
%------------------------------------------------------------------------
% Processing data (SCADA)
% Updating SCADA Voltage magnitudes: Basically freezing the current state
% SCADA Vm measurement until the next measurement arrives.
Number_of_timePoints_for_PMU_data    = PMU_SamplingFreq*SimulationTime;

Number_of_timePoints_for_SCADA_data  = length(SCADA_samples);

% Freezing SCADA until the next SCADA measurement is seen
N = Tar_SCADA_sr/Tar_PMU_sr;
SCADA_Vm_updated = zeros(size(PMU.Vm));
SCADA_Vm_updated(SCADA_samples/2,:) = SCADA.Vm;
for l = 1:Number_of_timePoints_for_SCADA_data-1
    SCADA_Vm_updated(SCADA_samples(l)/2:SCADA_samples(l+1)/2-1,:) = repmat(SCADA.Vm(l,:),N,1);
end
if l+1 == Number_of_timePoints_for_SCADA_data
    SCADA_Vm_updated(SCADA_samples(l+1)/2:Number_of_timePoints_for_PMU_data,:) = repmat(SCADA.Vm(l+1,:),N,1);
end
SCADA.Vm = SCADA_Vm_updated;

% Computing the SCADA Measurements [Real and Reactive powers]:
for i_bus = 1:n_bus
    
    if sum(SCADA_locations == i_bus)
        
        nodeBranches = [];
        nodeBranchIndx = [];
        for i_b = 1:length(AllBranches)
            if AllBranches(i_b,1) == i_bus || AllBranches(i_b,2) == i_bus
                nodeBranches = [nodeBranches; AllBranches(i_b,:)];
                nodeBranchIndx = [nodeBranchIndx; i_b];
            end
        end
        
        SCADA_P = zeros(Number_of_timePoints_for_SCADA_data,length(nodeBranchIndx));
        SCADA_Q = zeros(Number_of_timePoints_for_SCADA_data,length(nodeBranchIndx));
        SCADA_d = cell(1, length(nodeBranchIndx));
        
        for i_b = 1:size(nodeBranches,1)
            SCADA_P(:,i_b) = (real(fromBusPowInj(nodeBranchIndx(i_b),SCADA_samples)))';
            SCADA_Q(:,i_b) = (imag(fromBusPowInj(nodeBranchIndx(i_b),SCADA_samples)))';
            SCADA_d{i_b}   = sprintf('I_%d_%d',nodeBranches(i_b,1),nodeBranches(i_b,2));
        end
        
        SCADA.P{i_bus} = SCADA_P;
        SCADA.Q{i_bus} = SCADA_Q;
        SCADA.d{i_bus} = SCADA_d;
    else
        % If node doesn't have a SCADA sensor
        SCADA.P{i_bus} = [];
        SCADA.Q{i_bus} = [];
        SCADA.d{i_bus} = [];
    end
end

% Freezing SCADA P and Q variables until the next measurement arrives
n_SCADA = length(SCADA_locations);

for i_SCADA = 1:n_SCADA
    if isempty(SCADA.P{SCADA_locations(i_SCADA)})
        SCADA_P_updated = [];
        SCADA_Q_updated = [];
    else
        n_branch = size(SCADA.P{SCADA_locations(i_SCADA)},2);
        SCADA_P_updated = zeros(Number_of_timePoints_for_SCADA_data*N,n_branch);
        SCADA_Q_updated = zeros(Number_of_timePoints_for_SCADA_data*N,n_branch);
        for l = 1:Number_of_timePoints_for_SCADA_data
            SCADA_P_updated((l-1)*N+1:l*N,:) = repmat(SCADA.P{SCADA_locations(i_SCADA)}(l,:),N,1);
            SCADA_Q_updated((l-1)*N+1:l*N,:) = repmat(SCADA.Q{SCADA_locations(i_SCADA)}(l,:),N,1);
        end
    end
    SCADA.P{SCADA_locations(i_SCADA)} = SCADA_P_updated;
    SCADA.Q{SCADA_locations(i_SCADA)} = SCADA_Q_updated;
end

%------------------------------------------------------------------------
%{
% If the Euclidean norm based metrics needs to be computed using SCADA and
% PMU measurements [voltage magnitudes]
if compute_metrics
    % Taking the difference between the PMU measurements and SCADA measurements
    % to compute the norm based metrics
    % Initialization
    PMU_SCADA_difference = PMU.Vm(:,PMU_locations) - SCADA.Vm(:, PMU_locations);

    % Computing the norms 
    Types_of_metric = {1,2,'Inf'};
    PMU_SCADA_difference_norm = zeros(3,n_PMU);
    for i_PMU = 1:n_PMU
        for i_metric = 1:3
            PMU_SCADA_difference_norm(i_metric,i_PMU) = norm(PMU_SCADA_difference(:,i_PMU),Types_of_metric{i_metric});
        end
    end
    
    % Finding the metric data to be saved in CSV files
    metric_data = cell(2,3);
    for j_norm = 1:3
        [val, indx] = maxk(PMU_SCADA_difference_norm(j_norm,:),5);
        metric_data{1,j_norm} = num2str(PMU_locations(indx));
        metric_data{2,j_norm} = num2str(val(1));
    end

    scenDes_full{scenIdx_local+1,1}     = num2str(scenIdx);
    scenDes_full(scenIdx_local+1,2:14)  = scenDes(2,1:13);
    scenDes_full(scenIdx_local+1,15:18) = scenDes(2,14:17);
    scenDes_full(scenIdx_local+1,19:21) = metric_data(1,:);
    scenDes_full(scenIdx_local+1,22:24) = metric_data(2,:);
end 
%}
%------------------------------------------------------------------------
%% Plots
%------------------------------------------------------------------------
%% Real Power plots with error (PMU & SCADA)
%{
real_power_plot = figure;
subplot(3,1,1)
plot(t(PMU_samples),PMU.P{PMU_locations(2)},'LineWidth',2)
xlabel('time [s]')
title(sprintf('Real power at Bus %d [PMU]',PMU_locations(2)))
subplot(3,1,2)
plot(t(PMU_samples),SCADA.P{PMU_locations(2)},'LineWidth',2)
xlabel('time [s]')
title(sprintf('Real power at Bus %d [SCADA]',PMU_locations(2)))
subplot(3,1,3)
plot(t(PMU_samples),PMU.P{PMU_locations(2)}-SCADA.P{PMU_locations(2)},'LineWidth',2)
xlabel('time [s]')
%}
%------------------------------------------------------------------------
%% Voltage plots (PMU & SCADA)
%{
voltage_magnitude_plot = figure;
subplot(3,1,1)
plot(t(PMU_samples),PMU.Vm(:,PMU_locations), 'LineWidth',2)
xlabel('time [s]')
ylabel('Vm [pu]')
title('PMU measurements')
subplot(3,1,2)
plot(t(PMU_samples),SCADA.Vm(:,PMU_locations), 'LineWidth',2)
xlabel('time [s]')
ylabel('Vm [pu]')
title('SCADA measurements')
subplot(3,1,3)
plot(t(PMU_samples), PMU_SCADA_difference,'LineWidth',1)
xlabel('time [s]')
ylabel('$\Delta$Vm [pu]')
title('Difference between PMU and SCADA measurements')
%}
%------------------------------------------------------------------------
%% PMU Frequency plots 
frequency_plot = figure('Position',[-1680 528 560 420]);
freqs = PMU.f(:, PMU_locations);
mean_freq = mean(freqs,2);
plot(t(PMU_samples(1:end)), PMU.f(:, PMU_locations),'LineWidth',1.5); hold on
% plot(t(PMU_samples), 60.1*ones(size(mean_freq)), '--', 'color', color2,'LineWidth',1.5);
% plot(t(PMU_samples), 59.9*ones(size(mean_freq)), '--', 'color', color2,'LineWidth',1.5);
xlabel('time [s]');
ylabel('Frequency [Hz]');
title(sprintf('Scenario: %d', scenIdx));
% ylim([59, 61]);
%------------------------------------------------------------------------
%% PMU Voltage Magnitude and Angle plots
voltage_plot = figure('Position',[-1120 528 560 420]);
subplot(2,1,1)
plot(t(PMU_samples),PMU.Vm(:, PMU_locations), 'LineWidth',2)
title(sprintf('Scenario: %d', scenIdx));
xlabel('time [s]')
ylabel('Vm [pu]')
subplot(2,1,2)
plot(t(PMU_samples),PMU.Va(:, PMU_locations), 'LineWidth',2)
xlabel('time [s]');
ylabel('Voltage angle [deg]')
title(sprintf('Scenario: %d', scenIdx));
%}
%------------------------------------------------------------------------
%% Metric Plots
%{
impact_plot = figure;
subplot(3,1,1)
stem(PMU_SCADA_difference_norm(1,:),'LineWidth',1.5)
xlim([0, n_PMU+1])
xticks(1:n_PMU);
str_PMU_data = string(PMU_locations);
xticklabels(str_PMU_data);
xlabel('Location of PMU')
ylabel('Impact')
title('1-norm (average deviation)')
subplot(3,1,2)
stem(PMU_SCADA_difference_norm(2,:),'LineWidth',1.35)
xlim([0, n_PMU+1])
xticks(1:n_PMU);
xticklabels(str_PMU_data);
xlabel('Location of PMU')
ylabel('Impact')
title('2-norm (average deviation)')
subplot(3,1,3)
stem(PMU_SCADA_difference_norm(3,:),'LineWidth',1.2)
xlim([0, n_PMU+1])
xticks(1:n_PMU);
xticklabels(str_PMU_data);
xlabel('Location of PMU')
ylabel('Impact')
title('Infty-norm (worst-case deviation)');
%}
%------------------------------------------------------------------------
%% Real Power plots 
%{
if num_area == 1
    figure    
    plot(t, pelect, 'linewidth',1.5)
    xlabel('time [s]')
    title('Generator electrical powers [pu]')
elseif num_area == 2
    gen_buses = 53:68;
    area1_gens = gen_buses(area1_buses);
    area2_gens = gen_buses(area2_buses);
    figure
    subplot(2,1,1)
    plot(t, pelect(area1_buses,:), 'linewidth',1.5)
    xlabel('time [s]')
    title('Area-1 Generator electrical powers [pu]')
%     legend(num2str(area1_gens(:)),'Location','best','NumColumns',5)
%     legend boxoff
    subplot(2,1,2)
    plot(t, pelect(area2_buses,:), 'linewidth',1.5)
    xlabel('time [s]')
    title('Area-2 Generator electrical powers [pu]')
%     legend(num2str(area2_gens(:)),'Location','best','NumColumns',4)
%     legend boxoff
end
%}
%------------------------------------------------------------------------
%% AGC Plots 
if agc_control
    if num_area == 1
        ace_plot = figure;
        subplot(2,1,1)
%         stairs(t(ACE_data.k),ACE_data.ACE,'k-','linewidth',3);
        stairs(t(ACE_data.k),ACE_data.myACE,'k-','linewidth',3);
        xlabel('time [s]')
        title('Area Control Error')
        
        subplot(2,1,2)
        plot(t,tg_sig,'linewidth',2)
        xlabel('time [s]')
        title('Generator control inputs')
    elseif num_area == 2
        ace_plot = figure;
        subplot(3,1,1)
        stairs(t(ACE_data.k),ACE_data.myACE1,'k','linewidth',3,'displayname','ACE-1'); hold on
        stairs(t(ACE_data.k),ACE_data.myACE2, 'r','linewidth',3,'displayname','ACE-2')
        xlabel('time [s]')
        title('Area Control Error')
        legend location best
        legend box off
        
        subplot(3,1,2)
        plot(t,tg_sig(area1_buses,:),'linewidth',2)
        xlabel('time [s]')
        title('Area-1 generator control inputs')
        subplot(3,1,3)
        plot(t,tg_sig(area2_buses,:),'linewidth',2)
        title('Area-2 generator control inputs')
        xlabel('time [s]')
    end
end

% Down-sampling turbine governor control setpoints such that it matches
% with the PMU sampling frequency
% tg_sig = (tg_sig(:,1:2:end-1))'; % Rows - time; columns - generators 
%{
figure
plot(t, sum(tg_sig,1), 'linewidth',3)
xlabel('time [s]')
ylabel('sum of set points [pu]')
%}