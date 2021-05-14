function f = mtg_sig_PD(t,k)
% Syntax: f = mtg_sig(t,k)
% 12:37 PM 7/0/98
% defines modulation signal for turbine power reference
global tg_sig n_tg n_tgh
% fprintf('t = %0.2f;\t k = %d\n',t,k)
%==========================================================================
global busFreq TimeStep_of_simulation basrad
%==========================================================================
% Modified to incorporate AGC control in PowerDrone project
% Sai Pushpak Nandanoori - June 19, 2019
global tg_pot mac_spd
global agc_z agc_zm freq_zm pelect agc_ratio tg_con mac_con bus_int bus_v line
global agc_time_step agc_control agc_z1 agc_z2 agc_zm1 agc_zm2
global ace1 ace2 one_area S1_agcm S2_agcm num_area
global basmva TieLineScheduledPowers
global PMU_locations PMU_samples PMU
global PMU_area1_locations PMU_area2_locations
global ACE_data count area1_buses area2_buses
count = count + 1;

nominalFreq = basrad/2/pi;
B = 1.5929*10^4;
B2_area = B/2;
K = 1;
jay = sqrt(-1);
% The AGC implementation here is a simpler version of the AGC control
% usually done. 
% For one area AGC: tie-line power transfers [Not applicable]
% For two-area AGC: yet to model the tie-line power transfers in ACE
% computations. 
%==========================================================================
% for each iteration, this function is called three times.
f=0; %dummy variable
% fprintf('(t,k) = (%f,%d)\n',t,k);
% =========================================================================
if n_tg~=0||n_tgh~=0
    tg_sig(:,k) = zeros(n_tg+n_tgh,1);
    if agc_control       
%         disp('--------*--------')
%         if t>=30 && t <=43
%         fprintf('(t,k) = (%f,%d)\n',t,k);
%         fg = mod(t,agc_time_step);
%         fprintf('Resultant value = %0.4f\n',fg);
%         end
%         disp('--------*--------')
%         pause(0.4)
        if t<=0.0
            tg_sig(:,k) = zeros(n_tg+n_tgh,1);            
%         elseif mod(t,agc_time_step) == 0.0 
        elseif mod(k-1,agc_time_step*100) == 0.0 
            %disp('AGC loop')
            if num_area == 1               
                % ACE = mean(tg_con(:,3) - mac_spd(:,k)); %
                row_location = find(PMU_samples == k-1, 1);                
                myACE = mean(1-PMU.f(row_location,PMU_locations)/60);
                % ACE_data.ACE = [ACE_data.ACE, ACE];
                ACE_data.myACE = [ACE_data.myACE, myACE];
                ACE_data.k = [ACE_data.k, k];
                agc_z = agc_z+1*K*(-agc_z+B*myACE+sum(pelect(:,k)));
                agc_zm = [agc_zm agc_z];                
                tg_sig(:,k) = agc_z.*agc_ratio./mac_con(:,3)*basmva-tg_pot(:,5);                
            elseif num_area == 2 % two area                                
                % ACE1 = mean(tg_con(area1_buses,3) - mac_spd(area1_buses,k));
                % ACE2 = mean(tg_con(area2_buses,3) - mac_spd(area2_buses,k));
                row_location = find(PMU_samples == k-1, 1);
                %==========================================================================
                % Compute the tie line powers for ACE computation
                % Tie-Lines for IEEE 68 bus system
                % Line connecting buses 1-2; buses 1-27 and line 8-9
                % To compute the real power, we need the complex voltage
                % and complex current
%                 pol2cart(deg2rad(PMU.Va(:,i_bus)),PMU.Vm(:,i_bus));
                [V_real, V_imag] = pol2cart(deg2rad(PMU.Va(row_location,[1 1  9])),PMU.Vm(row_location,[1 1  9]));
                V1               = V_real' + jay*V_imag';
                [V_real, V_imag] = pol2cart(deg2rad(PMU.Va(row_location,[2 27 8])),PMU.Vm(row_location,[2 27 8]));
                V2               = V_real' + jay*V_imag';                
                R                = line([1 86 16],3); 
                X                = line([1 86 16],4); 
                B                = line([1 86 16],5);
                tap              = line([1 86 16],6); 
                phi              = line([1 86 16],7);
                
                [PQ_from_area1_to_area2, PQ_to]   = line_pq(V1, V2, R, X, B, tap, phi);
                [PQ_from_area2_to_area1, PQ_to]   = line_pq(V2, V1, R, X, B, tap, phi);
                %disp('--------*--------')
                %fprintf('(t,k) = (%f,%d)\n',t,k);
                % disp(PQ_from_area1_to_area2)
                % disp(PMU.f(row_location,:))
                %disp('--------*--------')
                tielinedifferences = (real(PQ_from_area1_to_area2) - TieLineScheduledPowers);
                ACE_data.tielinedifferences = [ACE_data.tielinedifferences, tielinedifferences];
                %{
                figure(15)
                subplot(2,1,1)
                plot(t,real(PQ_from_area1_to_area2),'o','MarkerSize',12); hold on                
                subplot(2,1,2)
                plot(t,PMU.f(row_location,:),'*','MarkerSize',12); hold on
                %}
                %{
                figure(16)
                plot(k, tielinedifferences, 'o'); hold on
                %}
                %==========================================================================
                % ACE computation:% Find the frequency and tie-line power errors:
                myACE1 = B2_area*(mean(1-(PMU.f(row_location, PMU_area1_locations)/60))) + 1*sum(real(PQ_from_area1_to_area2) - TieLineScheduledPowers);
                myACE2 = B2_area*(mean(1-(PMU.f(row_location, PMU_area2_locations)/60))) + 1*sum(real(PQ_from_area2_to_area1) + TieLineScheduledPowers);

%                 myACE1 = mean(1-(PMU.f(row_location, PMU_area1_locations)/60));
%                 myACE2 = mean(1-(PMU.f(row_location, PMU_area2_locations)/60));
                
                % ACE_data.ACE1 = [ACE_data.ACE1, ACE1];
                % ACE_data.ACE2 = [ACE_data.ACE2, ACE2];
                
                ACE_data.myACE1 = [ACE_data.myACE1, myACE1];
                ACE_data.myACE2 = [ACE_data.myACE2, myACE2];
                
                ACE_data.k = [ACE_data.k, k];                                
                agc_z1 = agc_z1+1*K*(-agc_z1+myACE1+sum(pelect(area1_buses,k)));
                agc_z2 = agc_z2+1*K*(-agc_z2+myACE2+sum(pelect(area2_buses,k)));                
                
                agc_zm1 = [agc_zm1 agc_z1];
                agc_zm2 = [agc_zm2 agc_z2];
                
                S1_agcm = agc_z1.*agc_ratio(area1_buses)./mac_con(area1_buses,3)*basmva-tg_pot(area1_buses,5);
                S2_agcm = agc_z2.*agc_ratio(area2_buses)./mac_con(area2_buses,3)*basmva-tg_pot(area2_buses,5);                
                tg_sig(:,k) = [S1_agcm; S2_agcm];
            end
        else
%             disp('--------*--------')
%             fprintf('(t,k) = (%f,%d)\n',t,k);
%             disp('--------*--------')
            tg_sig(:,k) = tg_sig(:,k-1);
        end
    end
end

return
