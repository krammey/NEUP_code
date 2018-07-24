clear all
close all
% Run this script first
start_time = tic;

%==========================================================================================
% Inputs
%==========================================================================================
% Salt Tank Inputs          |  Value            | Units       | Description/Notes
%---------------------------|-------------------|-------------|-----------------------------
TANK_RADIUS                 = 6;                % m             both tanks have same dimensions
TANK_HEIGHT                 = 6;                % m             in Milestone 3 report a tank height of 12 was used
STARTING_HOT_TANK_VOL_FRACT = 0.75;             % -             total fluid volume equals the volume of one tank.  This is the fraction of the total fluid that's initially in the hot tank
P_IN                        = 50E6;             % J/s           rate at which THERMAL power is transferred in when charging
P_OUT                       = 30E6;             % J/s           rate at which THERMAL power is transferred out when discharging
T_inf                       = 300;              % K             ambient temperature outside the tank
%------------------------------------------------------------------------------------------
% Sensitvity Analysis Input Ranges (9 combos used based on temp, hx_cf, and salt)
%------------------------------------------------------------------------------------------
T_HL                        = [550, 525, 500];  % K             array of target hot tank temps
T_CL                        = [350, 375, 400];  % K             array of target cold tank temps
hx_cf                       = [20, 35, 50];     %               heat transfer coeff for natural circulation convective heat xfer to environment - W/m2/K
salt                        = "Solar Salt";     %               Solar Salt, Hitec XL or Glauber
c_p                         = 0;                % J/kg-K        specific heat capacity , value initiated with 0 
rho                         = 0;                % kg/m3         density, value initiated with 0
if salt == "Solar Salt"
    c_p                     = 1496;
    rho                     = 1899;
elseif salt == "Hitec XL"
    c_p                     = 1447;
    rho                     = 1992;
elseif salt == "Glauber"
    c_p                     = 3349;
    rho                     = 1460;
else
    print("Error in salt input");
end
%------------------------------------------------------------------------------------------
% Sets up events (discharge, nothing, charge)
%------------------------------------------------------------------------------------------
N_EVENTS                    = 8;                %               number of events
event                       = [];               %               the event array is initiated with 0, 1 or 2 to correspond to a 
% charge, constant, or discharge cycle. can string together as many events 0, 1 and 2 as you want, 
% but code will give garbage results if either tank becomes totally empty or if heat loss to environment is 
% too high this can be fixed with some more coding but for now just keep an eye on the plots
endt_event                  = [];
tsteps_per_event            = [];
event(1)                    = 0;                % discharge
event(2)                    = 1;                % do nothing for a while
event(3)                    = 2;                % charge
event(4)                    = 1;
event(5)                    = 0;
event(6)                    = 2;
event(7)                    = 1;
event(8)                    = 0;
endt_event(1)               = 1800;             % s             event 1 ends after 1800 s (30 minutes)
endt_event(2)               = 3600;             % s             event 2 ends at 3600 s (hence lasts 1800 s)
endt_event(3)               = 5400;             % s             sets the end time of each event
endt_event(4)               = 12000;            % s
endt_event(5)               = 15000;            % s
endt_event(6)               = 17600;            % s
endt_event(7)               = 36000;            % s
endt_event(8)               = 38400;            % s
dt                          = 100;              % s             time step in seconds - I haven't tested convergence, could possibly jack this way up
%==========================================================================================
% end of inputs
%==========================================================================================




%-------------------------------------------------------------------------
% Initial Calculations
%-------------------------------------------------------------------------
TANK_TOTVOL             = pi*TANK_RADIUS*TANK_RADIUS*TANK_HEIGHT;       % m3, tank volume
TANK_SIDESA             = 2*pi*TANK_RADIUS*TANK_HEIGHT;                 % m2, area of tank sides
tsteps_per_event(1)     = round(endt_event(1)/dt);                      % number of time steps in the first event
endt_event(1)           = dt*tsteps_per_event(1);
nt                      = tsteps_per_event(1);
for i=2:N_EVENTS        % calculates total number of time steps for all events
    tsteps_per_event(i) = round((endt_event(i)-endt_event(i-1))/dt);% number of time steps in event i
    endt_event(i)       = endt_event(i-1) + dt*tsteps_per_event(i); % rounds entries here?
    nt                  = nt + tsteps_per_event(i);
end
t                       = (0:dt:(nt-1)*dt); % time array?

%-------------------------------------------------------------------------
% Creates arrays of length nt for time, temperatures, volumes, and masses
%-------------------------------------------------------------------------
Th                      = zeros(1,nt);
Tc                      = zeros(1,nt);
Vh                      = zeros(1,nt);
Vc                      = zeros(1,nt);
mh                      = zeros(1,nt);
mc                      = zeros(1,nt);
% mdotchg = zeros(nt,1)';
% mdotdis = zeros(nt,1)';

% Creates temporary 2D arrays for Th, Tc, Vh, and Vc
% These will store Th, Tc, Vh, and Vc for all nine sensitivity trials  
Th_temp                 = zeros(9,nt);
Tc_temp                 = zeros(9,nt);
Vh_temp                 = zeros(9,nt);
Vc_temp                 = zeros(9,nt);



% These two loops will iterate through the different sensitivity inputs
count = 1;
for h=1:3       % h iterates through the heat transfer coefficients
    for g=1:3   % g iterates though the temperature ranges 
        % inital assignments/calculations for Th, Tc, Vc, and Vh at time step zero
        Th(1) = T_HL(g); % initial T in hot tank
        Tc(1) = T_CL(g); % initial T in cold tank
        mh(1) = TANK_TOTVOL*rho*STARTING_HOT_TANK_VOL_FRACT;     % initial mass in hot tank
        mc(1) = TANK_TOTVOL*rho*(1-STARTING_HOT_TANK_VOL_FRACT); % initial mass in cold tank
        Vh(1) = mh(1)/rho; % initial volume in hot tank
        Vc(1) = mc(1)/rho; % initial volume in cold tank
        % mdotchg(1) = 0;
        % mdotdis(1) = 0;
        ix = 2;
        mt = [];
        % calculated temperatures and volumes at the time steps for each case
        % the formulas used for Th, Tc, Vh, and Vc can be found in the report -- derived by Schneider 
        for i=1:N_EVENTS
            if event(i)==0  % discharge
                for j=1:tsteps_per_event(i) 
                    Th(ix) = Th(ix-1)-dt*(hx_cf(h)*(Th(ix-1)-T_inf)*TANK_SIDESA/(c_p*rho*Vh(ix-1)));
                    Tc(ix) = Tc(ix-1)-dt*(hx_cf(h)*(Tc(ix-1)-T_inf)*TANK_SIDESA/(c_p*rho*Vc(ix-1)));
                    Vh(ix) = Vh(ix-1)-dt*(P_OUT/(rho*c_p*(Th(ix-1)-T_CL(g)))); 
                    Vc(ix) = Vc(ix-1)+dt*(P_OUT/(rho*c_p*(Th(ix-1)-T_CL(g))));
                    mt     = dt*(P_OUT/(rho*c_p*(Th(ix-1)-T_CL(g)))); % Why is this the same as the second term in V eqn?********
                    ix     = ix+1;
                end % for int j
                ix = ix-2; % reset time step index after each event
            elseif event(i) == 1  % just heat loss
                for j=1:tsteps_per_event(i)
                     Th(ix) = Th(ix-1)-dt*(hx_cf(h)*(Th(ix-1)-T_inf)*TANK_SIDESA/(c_p*rho*Vh(ix-1)));
                     Tc(ix) = Tc(ix-1)-dt*(hx_cf(h)*(Tc(ix-1)-T_inf)*TANK_SIDESA/(c_p*rho*Vc(ix-1)));
                     Vh(ix) = Vh(ix-1);
                     Vc(ix) = Vc(ix-1);
                     ix     = ix+1;
                end %for int j
                ix = ix-2; % reset time step index after each event
            elseif event(i) == 2  % charge
                for j=1:tsteps_per_event(i)
                    Th(ix) = Th(ix-1)+dt*(P_IN*(T_HL(g)-Th(ix-1))/(Vh(ix-1)*c_p*rho*(T_HL(g)-Tc(ix-1)))-hx_cf(h)*(Th(ix-1)-T_inf)*TANK_SIDESA/(c_p*rho*Vh(ix-1)));
                    Tc(ix) = Tc(ix-1)-dt*(P_IN*(T_HL(g)-Th(ix-1))/(Vh(ix-1)*c_p*rho*(T_HL(g)-Tc(ix-1)))+hx_cf(h)*(Tc(ix-1)-T_inf)*TANK_SIDESA/(c_p*rho*Vc(ix-1)));
                    Vh(ix) = Vh(ix-1)+dt*(P_IN/(rho*c_p*(T_HL(g)-Tc(ix-1))));
                    Vc(ix) = Vc(ix-1)-dt*(P_IN/(rho*c_p*(T_HL(g)-Tc(ix-1)))); 
                    ix     = ix+1;
                end %for int j
                ix = ix-2; % reset time step index after each event
            else
                break;
            end
        end % for int i
        % once all temperatures and volumes are calculated for a given sensitivity trial (1-9), the results are stored
        % in the 2D temporary arrays
        for r = 1:nt
            Th_temp(count,r) = Th(r);
            Tc_temp(count,r) = Tc(r);
            Vh_temp(count,r) = Vh(r);
            Vc_temp(count,r) = Vc(r);
        end %for int r

        count = count+1;
    
    end % for int g
    
end %for int h

%these 1D arrays are created to hold the temperature and volume results of each sensitivity trial
%9 possible combos for Th, Tc, Vh, and Vc
Th1 = zeros(1,nt);
Tc1 = zeros(1,nt);
Vh1 = zeros(1,nt);
Vc1 = zeros(1,nt);

Th2 = zeros(1,nt);
Tc2 = zeros(1,nt);
Vh2 = zeros(1,nt);
Vc2 = zeros(1,nt);

Th3 = zeros(1,nt);
Tc3 = zeros(1,nt);
Vh3 = zeros(1,nt);
Vc3 = zeros(1,nt);

Th4 = zeros(1,nt);
Tc4 = zeros(1,nt);
Vh4 = zeros(1,nt);
Vc4 = zeros(1,nt);

Th5 = zeros(1,nt);
Tc5 = zeros(1,nt);
Vh5 = zeros(1,nt);
Vc5 = zeros(1,nt);

Th6 = zeros(1,nt);
Tc6 = zeros(1,nt);
Vh6 = zeros(1,nt);
Vc6 = zeros(1,nt);

Th7 = zeros(1,nt);
Tc7 = zeros(1,nt);
Vh7 = zeros(1,nt);
Vc7 = zeros(1,nt);

Th8 = zeros(1,nt);
Tc8 = zeros(1,nt);
Vh8 = zeros(1,nt);
Vc8 = zeros(1,nt);

Th9 = zeros(1,nt);
Tc9 = zeros(1,nt);
Vh9 = zeros(1,nt);
Vc9 = zeros(1,nt);

% assignments from the 2D temporary array are made to the individual 1D arrays defined above 
for y=1:nt
    Th1(y) = Th_temp(1,y);
    Tc1(y) = Tc_temp(1,y);
    Vh1(y) = Vh_temp(1,y);
    Vc1(y) = Vc_temp(1,y);

    Th2(y) = Th_temp(2,y);
    Tc2(y) = Tc_temp(2,y);
    Vh2(y) = Vh_temp(2,y);
    Vc2(y) = Vc_temp(2,y);

    Th3(y) = Th_temp(3,y);
    Tc3(y) = Tc_temp(3,y);
    Vh3(y) = Vh_temp(3,y);
    Vc3(y) = Vc_temp(3,y);

    Th4(y) = Th_temp(4,y);
    Tc4(y) = Tc_temp(4,y);
    Vh4(y) = Vh_temp(4,y);
    Vc4(y) = Vc_temp(4,y);

    Th5(y) = Th_temp(5,y);
    Tc5(y) = Tc_temp(5,y);
    Vh5(y) = Vh_temp(5,y);
    Vc5(y) = Vc_temp(5,y);

    Th6(y) = Th_temp(6,y);
    Tc6(y) = Tc_temp(6,y);
    Vh6(y) = Vh_temp(6,y);
    Vc6(y) = Vc_temp(6,y);

    Th7(y) = Th_temp(7,y);
    Tc7(y) = Tc_temp(7,y);
    Vh7(y) = Vh_temp(7,y);
    Vc7(y) = Vc_temp(7,y);

    Th8(y) = Th_temp(8,y);
    Tc8(y) = Tc_temp(8,y);
    Vh8(y) = Vh_temp(8,y);
    Vc8(y) = Vc_temp(8,y);

    Th9(y) = Th_temp(9,y);
    Tc9(y) = Tc_temp(9,y);
    Vh9(y) = Vh_temp(9,y);
    Vc9(y) = Vc_temp(9,y);
end % for y


%==========================================================================================
% Plotting
%==========================================================================================

%------------------------------------------------------------------------------------------
% Temperature Plot
%------------------------------------------------------------------------------------------
subplot(1,2,1)
plot(t,Th9,t,Th8,t,Th7,t,Th6,t,Th5,t,Th4,t,Th3,t,Th2,t,Th1,...
    t,Tc9,'-.',t,Tc8,'-.',t,Tc7,'-.',t,Tc6,'-.',t,Tc5,'-.',t,Tc4,'-.',t,Tc3,'-.',t,Tc2,'-.',t,Tc1,'-.')
legend('9','8','7','6','5','4','3','2','1','9','8','7','6','5','4','3','2','1')
xlim([0 max(t)]), xlabel('Time (s)'), ylabel('Temperature [K]')
set(gcf,'Color','w')

%------------------------------------------------------------------------------------------
% Volume Plot
%------------------------------------------------------------------------------------------
subplot(1,2,2)
plot(t,Vh9,t,Vh8,t,Vh7,t,Vh6,t,Vh5,t,Vh4,t,Vh3,t,Vh2,t,Vh1,...
    t,Vc9,'-.',t,Vc8,'-.',t,Vc7,'-.',t,Vc6,'-.',t,Vc5,'-.',t,Vc4,'-.',t,Vc3,'-.',t,Vc2,'-.',t,Vc1,'-.')
legend('9','8','7','6','5','4','3','2','1','9','8','7','6','5','4','3','2','1')
xlim([0 max(t)]), xlabel('Time (s)'), ylabel('Volume [m^3]')
set(gcf,'Color','w')


%------------------------------------------------------------------------------------------
% All Separate Plots
%------------------------------------------------------------------------------------------
figure
subplot(2,2,1), title('Hot Tank Temperature')
plot(t,Th9,t,Th8,t,Th7,t,Th6,t,Th5,t,Th4,t,Th3,t,Th2,t,Th1)
legend('9','8','7','6','5','4','3','2','1')
xlim([0 max(t)]), xlabel('Time (s)'), ylabel('Temperature [K]')
set(gcf,'Color','w')

subplot(2,2,2), title('Hot Tank Volume')
plot(t,Vh9,t,Vh8,t,Vh7,t,Vh6,t,Vh5,t,Vh4,t,Vh3,t,Vh2,t,Vh1)
legend('9','8','7','6','5','4','3','2','1')
xlim([0 max(t)]), xlabel('Time (s)'), ylabel('Volume [m^3]')
set(gcf,'Color','w')

subplot(2,2,3), title('Cold Tank Temperature')
plot(t,Tc9,'-.',t,Tc8,'-.',t,Tc7,'-.',t,Tc6,'-.',t,Tc5,'-.',t,Tc4,'-.',t,Tc3,'-.',t,Tc2,'-.',t,Tc1,'-.')
legend('9','8','7','6','5','4','3','2','1')
xlim([0 max(t)]), xlabel('Time (s)'), ylabel('Temperature [K]')
set(gcf,'Color','w')

subplot(2,2,4), title('Cold Tank Volume')
plot(t,Vc9,'-.',t,Vc8,'-.',t,Vc7,'-.',t,Vc6,'-.',t,Vc5,'-.',t,Vc4,'-.',t,Vc3,'-.',t,Vc2,'-.',t,Vc1,'-.')
legend('9','8','7','6','5','4','3','2','1')
xlim([0 max(t)]), xlabel('Time (s)'), ylabel('Volume [m^3]')
set(gcf,'Color','w')


% figure
% plot(t,Th_temp',t,Tc_temp','-.')



