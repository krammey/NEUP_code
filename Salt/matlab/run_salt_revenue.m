
caseNumber = 3;
%===========================================================================================
% Salt Tank properties
%===========================================================================================
% Parameter      |  Value                       | Units       | Description
%----------------|------------------------------|-------------|-----------------------------
D_RAMP_RATE         = 1.67;                     % percent/min,  discharge ramp rate (SA turbine, % of max power/min)
POWER_INITIAL       = 0;                        % MW,           cold start
DT                  = 10;  						% s,            time step

%===========================================================================================
% Main plant properties
%===========================================================================================
% Parameter      |  Value                       | Units       | Description
%----------------|------------------------------|-------------|-----------------------------
MAIN_POWER          = 1300;                     % MWe,          power from main turbine of plant
MIN_LEVEL           = 25;                       % percent,      minimum turbine level as percent of main power
MIN_LOAD            = (MIN_LEVEL/100)*MAIN_POWER;% MWe

%===========================================================================================
% Sinusoidal price curve, amortization values, and other economic info
%===========================================================================================
% Parameter      |  Value         |  Units       |  Description
%----------------|----------------|--------------|------------------------------------------
life                = 40;         % years,          amortization period
interest            = 0.07;       %                 for amortization period
period              = 6;          % hours,          price period
peakAmplitude       = 25;         % $/MWh
avgElecPrice        = 34;         % $/MWh
coldCyclesPerYear   = 100;        % cycles/year
warmCyclesPerYear   = 50;         % cycles/year
hotCyclesPerYear    = 50;         % cycles/year
var_om              = 8;          % $/MWh,          from Neal

% Define power and storage time inputs
POWER_array         = [50 100 150]; % array of powers to be tested
T_MAX_array         = [1 2 3 4 5 6 7 8]; % array of storage times to be tested

% Initialize results arrays and counter
Revenue_Results     = [POWER_array' zeros(length(POWER_array),length(T_MAX_array))];
annualCC_Results    = [POWER_array' zeros(length(POWER_array),length(T_MAX_array))];
TotalCC_Results     = [POWER_array' zeros(length(POWER_array),length(T_MAX_array))];
TotalOM_Results     = [POWER_array' zeros(length(POWER_array),length(T_MAX_array))];

count = 0;
for i = 1:length(POWER_array)
    hr_count = 0;
    for j = 1:length(T_MAX_array)
        hr_count = hr_count+1;
        count = count+1;
        POWER = POWER_array(i);
        T_MAX = T_MAX_array(j);
        
        %===========================================================================================
        % Power- and energy-dependent steam accumulator properties
        %===========================================================================================
        ENERGY   = POWER*T_MAX;
        T_RAMP   = 60*((POWER-POWER_INITIAL)/((D_RAMP_RATE/100)*POWER)); % s, time to ramp up to max power
        T_END    = T_RAMP+T_MAX*3600;                 % [s]  discharge time with ramp up and time at max power
        t_length = floor(T_RAMP/DT);                  %      length of ramping power vector
        t        = 1:t_length;                        %      time array
        pow      = DT*t*((D_RAMP_RATE/100)*POWER)/60; % [MW] power as SA turbine is ramping up
        
        %------------------------------------------------------------------------
        % Run revenue_model.m
        %------------------------------------------------------------------------
        addpath('../../Revenue/')
        disp(['Power = ' num2str(POWER) ' Hours = ' num2str(T_MAX)])
        [netRevenue,annualCC,startCost,annualOM,totalCC,fixedOM] = revenue_model("salt",[],T_END,POWER,ENERGY,...
            MAIN_POWER,MIN_LOAD,life,interest,period,peakAmplitude,...
            avgElecPrice,caseNumber,hotCyclesPerYear,warmCyclesPerYear,coldCyclesPerYear,var_om);
        disp(' ')
        
        %         %------------------------------------------------------------------------
        %         % Save Results
        %         %------------------------------------------------------------------------
        %         Revenue_Results(i,j+1) = netRevenue;
        %         CC_Results(i,j+1)      = CC;
        %         TotalCC_Results(i,j+1) = totalCC;
        %         TotalOM_Results(i,j+1) = totalOM;
        
        %------------------------------------------------------------------------
        % Save Results and convert units
        %------------------------------------------------------------------------
        Revenue_Results(i,j+1)  = netRevenue*1000/(ENERGY);
        annualCC_Results(i,j+1) = annualCC*1000/(ENERGY); % $/kWh(e)
        TotalCC_Results(i,j+1)  = totalCC*1000/(ENERGY);  % $/kWh
        TotalOM_Results(i,j+1)  = annualOM*1000/(ENERGY); % $/kWh(e)
        fixedOM_Results(i,j+1)  = fixedOM; % $/kW-year
    end
end



% ------------------------------------------------------------------------
% Plotting
% ------------------------------------------------------------------------
addpath('/Users/trins/Documents/MATLAB/mtit')
addpath('/Users/trins/Documents/MATLAB/altmany-export_fig-9d97e2c')
% COLORS
% Blue: 0 0.45 0.74
% grey = '0.4 0.4 0.4';
% Red: 0.85 0.33 0.1
red = "0.85 0.33 0.1";
blue = "0 0.45 0.74";
grey = "0.4 0.4 0.4";
MARKERS = ["-d";"-s";"-o"];
COLORS = [blue; red; grey];
% LABELS = ["Annual Net Revenue [million $]"; "Annual CC [million $]"; "Total O&M [million $]"; "Total CC [million $]"];
LABELS = ["Annual Net Revenue [$/kWh_e]"; "Annual CC [$/kWh_e]"; "Total O&M [$/kWh_e]"; "Total CC [$/kWh_e]";"Fixed O&M [$/kW-year]"]; % new labels after converting units

% Re-organize results to be able to plot using for loop 
Results        = zeros(size(Revenue_Results,1),size(Revenue_Results,2),0);
Results(:,:,1) = Revenue_Results;
Results(:,:,2) = annualCC_Results;
Results(:,:,3) = TotalOM_Results;
Results(:,:,4) = TotalCC_Results;
Results(:,:,5) = fixedOM_Results;

figure
for r=1:5 % loop through results: 1-Revenue, 2-CC, 3-Total OM, 4-Total CC
    subplot(2,3,r)
    df = Results(:,:,r);
    hold on
    for idx = 1:length(POWER_array)
        plot(categorical(T_MAX_array),df(idx,2:length(T_MAX_array)+1),MARKERS(idx),'MarkerFaceColor',COLORS(idx),'Color',COLORS(idx))
    end
    ylabel(LABELS(r)), xlabel('Hours of storage')
    hold off
end
legend(string(POWER_array),'Location','BestOutside'), legend('boxoff'),
set(gcf,'Color','w')
set(findall(gcf,'-property','FontSize'),'FontSize',10)
set(findall(gcf,'-property','MarkerSize'),'MarkerSize',4)
mtit(['Case Number = ' num2str(caseNumber)])

% export_fig(['Case' num2str(caseNumber) '.fig'])
disp(' ')

