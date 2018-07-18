
function [netRevenue,CC,RC,RD,totalOM,totalCC] = steam_revenue(object,POWER_SA,ENERGY_SA,MAIN_POWER,MIN_LOAD,~,life,interest,period,peakAmplitude,avgElecPrice,caseNumber, hotCyclesPerYear, warmCyclesPerYear, coldCyclesPerYear, var_om)   
    % length must be connected to capital cost through pipe and insulation costs'
    caseName = 'CASE' + string(caseNumber);
    Y   = 24*365/period;                % storage cycles per year
    d_t = object.discharge_time/3600;   % hours, discharge time
    c_t = object.charge_time/3600;      % hours, charge time for SA
    
    % Read in data with rows:
    % 1 caseNumber
    % 2 Pref
    % 3 Eref
    % 4 cc_p
    % 5 cc_e
    % 6 om_p
    % 7 om_e
    % 8 cc_p_scale
    % 9 cc_e_scale 
    % 10 om_p_scale 
    % 11 om_e_scale 
    % 12 coldStart 
    % 13 warmStart 
    % 14 hotStart
    data = xlsread('steam_revenue_input.xls');
    column = caseNumber+1;
    
    %reference power and energy assigned
    Pref = data(2,column); % reference power
    Eref = data(3,column); % reference energy
    cc_p = data(4,column); % cost (power)
    cc_e = data(5,column); % cost (energy)
    om_p = data(6,column); % OM (power)
    om_e = data(7,column); % OM (energy)

    %scale factors assigned for power/energy cost and OM 
    cc_p_scale = data(8,column); % scaleFactor cost (power)
    cc_e_scale = data(9,column); % scaleFactor cost (energy)
    om_p_scale = data(10,column); % scaleFactor OM (power)
    om_e_scale = data(11,column);% scaleFactor cost (energy)

    %Cp Ce Op Oe scaled and assinged 
    Cp_scaled = cc_p * (POWER_SA/Pref)^cc_p_scale;  %cost scaled (power)
    Ce_scaled = cc_e * (ENERGY_SA/Eref)^cc_e_scale; %cost scaled (energy)
    Op_scaled = om_p * (POWER_SA/Pref)^om_p_scale;  %OM scaled (power)
    Oe_scaled = om_e * (ENERGY_SA/Eref)^om_e_scale; %OM scaled (energy)

    %cycles OM calcaulated 
    coldStart   = data(12,column); %$/MW-Cycle
    warmStart   = data(13,column); %$/MW-Cycle
    hotStart    = data(14,column); %$/MW-Cycle
    cyclingCost = ((coldStart * coldCyclesPerYear + warmStart * warmCyclesPerYear + hotStart * hotCyclesPerYear) * POWER_SA)/1000000; %MM$/year

    %calculations to determine ADC, ACP, DP
    c1   = (3/4)*period-c_t/2;  % hr, charge time integral lower bound
    c2   = (3/4)*period+c_t/2;  % hr, charge time integral upper bound
    d1   = (period/4)-d_t/2;    % hr, discharge time integral lower bound
    d2   = (period/4)+d_t/2;    % hr, discharge time integral upper bound
    y    = @(t)peakAmplitude*sin((2*pi()*t)/period)+avgElecPrice; 
    intC = integral(y,c1,c2);   % $/MW
    intD = integral(y,d1,d2);   % $/MW
    ADP  = intD/(d2-d1);        % $/MWh, Average discharge price
    ACP  = intC/(c2-c1);        % $/MWh, Average charge price
    %DP   = ADP-ACP ;            % $/MWh, delta price

    % Final costs/revenues caluclated and displayed (These need to be updated to account for ramping.)
    totalCC    = Cp_scaled+Ce_scaled;
    disp(['Total overnight CC is ' num2str(totalCC) ' MM$'])
    fixed_om   = (Oe_scaled+Op_scaled)*(10^6)*(1/POWER_SA)*(1/1000);
    disp(['Fixed O&M is ' num2str(fixed_om) ' $/kw-year'])
    totalOM    = Op_scaled + Oe_scaled + cyclingCost + var_om *ENERGY_SA*Y/1000000;
    disp(['Total O&M is ' num2str(totalOM) ' MM$/year'])
    startCost  = cyclingCost*1000000/(POWER_SA*(hotCyclesPerYear+warmCyclesPerYear+coldCyclesPerYear));
    disp(['Start cost is ' num2str(startCost) ' $/MW-start'])
    RC         = ACP*c_t*Y*(MAIN_POWER-MIN_LOAD)/10^6;                % forgone revenue from charging
    RD         = ADP*d_t*Y*POWER_SA/10^6;                             % forgone revenue from discharging
    CC         = totalCC*(interest+(interest/((1+interest)^life-1))); % amortized capital cost
    netRevenue = RD-RC-CC-totalOM; % revenue provided by the addition of the accumulator
    disp(['Net revenue is ' num2str(netRevenue) ' MM$/year'])
    
%     disp(["Total overnight CC is 95.2914 MM$" ...
%     "Fixed O&M is 95.351 $/kw-year" ... %****
%     "Total O&M is 10.0495 MM$/year" ... %****
%     "Start cost is 61 $/MW-start" ...
%     "Net revenue is -582.1025 MM$/year" ... %****
%     "Total run time = 5.12 seconds."]')
end

