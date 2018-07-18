
function [netRevenue,CC,RC,RD,totalOM,totalCC] = revenue_model(TECH,object,T_DISCHARGE,POWER,ENERGY,MAIN_POWER,MIN_LOAD,~,life,interest,period,peakAmplitude,avgElecPrice,caseNumber, hotCyclesPerYear, warmCyclesPerYear, coldCyclesPerYear, var_om)   
    
    
    caseName = 'CASE' + string(caseNumber);
    Y   = 24*365/period;                           % storage cycles per year
    if TECH == "steam"
        FileName = 'steam_revenue_input.xls';      % file with inputs
        d_t = T_DISCHARGE/3600;                    % hrs, discharge time
        c_t = object.charge_time/3600;             % hrs, charge time for steam
    elseif TECH == "salt"
        FileName = 'salt_revenue_input.xls';
        d_t = T_DISCHARGE/3600;                    % hrs, discharge time
        c_t = (POWER+7.61232)/(POWER-7.61232)*d_t; % hrs, charge time for salt
    else 
        disp('ERROR: Set TECH input to either "steam" or "salt"')
        return
    end
    
    % Read in data with rows:
    % (1) caseNumber            (8) cc_p_scale
    % (2) Pref                  (9) cc_e_scale 
    % (3) Eref                  (10) om_p_scale 
    % (4) cc_p                  (11) om_e_scale 
    % (5) cc_e                  (12) coldStart 
    % (6) om_p                  (13) warmStart 
    % (7) om_e                  (14) hotStart
    data = xlsread(FileName);
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
    Cp_scaled = cc_p * (POWER/Pref)^cc_p_scale;  %cost scaled (power)
    Ce_scaled = cc_e * (ENERGY/Eref)^cc_e_scale; %cost scaled (energy)
    Op_scaled = om_p * (POWER/Pref)^om_p_scale;  %OM scaled (power)
    Oe_scaled = om_e * (ENERGY/Eref)^om_e_scale; %OM scaled (energy)

    %cycles OM calcaulated 
    coldStart   = data(12,column); %$/MW-Cycle
    warmStart   = data(13,column); %$/MW-Cycle
    hotStart    = data(14,column); %$/MW-Cycle
    cyclingCost = ((coldStart * coldCyclesPerYear + warmStart * warmCyclesPerYear + hotStart * hotCyclesPerYear) * POWER)/1000000; %MM$/year

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
    fixed_om   = (Oe_scaled+Op_scaled)*(10^6)*(1/POWER)*(1/1000);
    disp(['Fixed O&M is ' num2str(fixed_om) ' $/kw-year'])
    totalOM    = Op_scaled + Oe_scaled + cyclingCost + var_om *ENERGY*Y/1000000;
    disp(['Total O&M is ' num2str(totalOM) ' MM$/year'])
    startCost  = cyclingCost*1000000/(POWER*(hotCyclesPerYear+warmCyclesPerYear+coldCyclesPerYear));
    disp(['Start cost is ' num2str(startCost) ' $/MW-start'])
    RC         = ACP*c_t*Y*(MAIN_POWER-MIN_LOAD)/10^6;                % forgone revenue from charging
    RD         = ADP*d_t*Y*POWER/10^6;                             % forgone revenue from discharging
    CC         = totalCC*(interest+(interest/((1+interest)^life-1))); % amortized capital cost
    netRevenue = RD-RC-CC-totalOM; % revenue provided by the addition of the accumulator
    disp(['Net revenue is ' num2str(netRevenue) ' MM$/year'])

end

