
% Run this script first
% To test against old code versions run:
% run ./v2/run_steam_model
% run ./v1/RUN_SA_SEPARATE.m
% run ./v0/RUN_SA_SEPARATE.m
%===========================================================================================
clear all
close all
start_time = tic;
%===========================================================================================
% Steam Accumulator properties
%===========================================================================================
% Parameter      |  Value                       | Units       | Description
%----------------|------------------------------|-------------|-----------------------------
POWER               = 50;                       % MW,           maximum electric power provided by steam accumulator
ENERGY              = 400;                      % MWh 
T_MAX               = ENERGY/POWER;             % hr,           time capacity at max SA power --> essentially the discharge time without ramping
D_RAMP_RATE         = 1.67;                     % percent/min,  discharge ramp rate (SA turbine, % of max power/min)
POWER_INITIAL       = 0;                        % MW,           cold start
T_RAMP              = 60*((POWER-POWER_INITIAL)/((D_RAMP_RATE/100)*POWER)); % s, time to ramp up to max power
T_END               = T_RAMP+T_MAX*3600; 		% s,            discharge time with ramp up and time at max power
T_STORE             = 10*3600;                  % s,            storage time (10 hours)
DT                  = 10;  						% s,            time step
RTANK               = 0.4064;                   % m,            pipe radius (16 inches)
LTANK               = 100000.;                  % m,            pipe length
VTANK               = LTANK.*pi.*RTANK^2.;  	% m3
t_length            = floor(T_RAMP/DT);         %               length of ramping power vector
pow                 = zeros(t_length,1);		% MW,           power as SA turbine is ramping up
for t = 1:t_length
    pow(t)          = t*DT*((D_RAMP_RATE/100)*POWER)/60;
end

%===========================================================================================
% Accumulator initial thermo properties
%===========================================================================================
P0                  = 70; 						% bar,          initial pressure
X0                  = 0.06;                     %               vapor quality (mass fraction) 

%===========================================================================================
% Main plant properties
%===========================================================================================
% Parameter      |  Value                       | Units       | Description
%----------------|------------------------------|-------------|-----------------------------
MAIN_POWER          = 1300;                     % MWe,          power from main turbine of plant
THERMAL_POWER       = 3500;                     % MWt,          from steam generator
MDOTBASE            = 1333;                     % kg/s,         mass flow rate at main_power
MIN_LEVEL           = 25;                       % percent,      minimum turbine level as percent of main power
MIN_LOAD            = (MIN_LEVEL/100)*MAIN_POWER;% MWe
p_topup             = 1;                        % bar,          low pressure makeup tank
h_topup             = 206;                      % kJ/kg,        low pressure makeup tank
sgh_output          = 2770;                     % kJ/kg,        outlet enthalpy at steam gen
sgh_input           = sgh_output-(THERMAL_POWER*1000)/MDOTBASE; % kJ/kg, inlet enthalpy at steam gen
Number              = 1;                        %               number of cycles
POWER_REDUCTION     = MAIN_POWER-MIN_LOAD;
MDOT_CHARGE         = MDOTBASE-(MIN_LEVEL/100)*MDOTBASE;% kg/s, maximum charging mass flow rate to produce min_load at base case values

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

%------------------------------------------------------------------------
% Set Case Number
%------------------------------------------------------------------------
% DESIGNS: 
% divert main steam (MS)
% preheat feedwater (FW)
% POWER TRAIN: 
% if Pdisch > 0.2Preactor, need additional power train (PT) 
% if not, can just upgrade exisitng (UG)
% HEAT SINK: 
% if diverting MS and Pchg > 0.5Preactor, need heat sink in 
% case of issue taking accumulator offline (HS); 
% otherwise (NA)
% case 1: MS, PT, HS
% case 2: MS, UG, HS
% case 3: MS, PT, NA
% case 4: MS, UG, NA
% case 5: FW, PT, NA
% case 6: FW, PT, NA
caseNumber = 3;

%------------------------------------------------------------------------
% Run steam_model.m (onlyun either block 1 or block 2)
%------------------------------------------------------------------------
% Block 1 - Discharge cycle
discharge = 1; 										% turns discharge on
ACC = steam_model(T_END,T_STORE,T_RAMP,DT,VTANK,P0,X0,LTANK,RTANK,discharge);
ACC.run_accumulator(POWER, pow, discharge); 		% runs discharge
ACC.charge(P0,X0,VTANK,MDOT_CHARGE,POWER_REDUCTION,POWER,p_topup,h_topup,sgh_input,sgh_output);

% Block 2 - storage mode
% discharge = 0; % turns discharge off (store)
% acc=steam_model(T_END,T_STORE,T_RAMP,DT,VTANK,P0,X0,LTANK,RTANK,discharge);
% acc.run_accumulator(POWER_SA, pow, discharge);    % store

%------------------------------------------------------------------------
% Run revenue_model.m
%------------------------------------------------------------------------
addpath('../../Revenue/')
[netRevenue,CC,RC,RD,totalOM,totalCC] = revenue_model("steam",ACC,T_END,POWER,ENERGY,...
    MAIN_POWER,MIN_LOAD,LTANK,life,interest,period,peakAmplitude,...
    avgElecPrice,caseNumber,hotCyclesPerYear,warmCyclesPerYear,coldCyclesPerYear,var_om);
        
%------------------------------------------------------------------------
% Plotting  (need to add more plots)
%------------------------------------------------------------------------
% acc.plots();

fprintf('Total run time = %.2f seconds.\n', toc(start_time));
disp(' ')


