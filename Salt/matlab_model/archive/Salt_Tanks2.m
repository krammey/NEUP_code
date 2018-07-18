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
        
        % once all temperatures and volumes are calculated for a given trial, the results are stored in 2D temporary arrays
        Th_temp(count,:) = Th;
        Tc_temp(count,:) = Tc;
        Vh_temp(count,:) = Vh;
        Vc_temp(count,:) = Vc;
        count = count+1;
        
    end % for int g
    
end %for int h

%these 1D arrays are created to hold the temperature and volume results of each sensitivity trial
%9 possible combos for Th, Tc, Vh, and Vc


%==========================================================================================
% Plotting
%==========================================================================================

%------------------------------------------------------------------------------------------
% Temperature Plot
%------------------------------------------------------------------------------------------
figure
subplot(1,2,1)
plot(t,Th_temp',t,Tc_temp','-.')
legend('1','2','3','4','5','6','7','8','9','1','2','3','4','5','6','7','8','9')
xlim([0 max(t)]), xlabel('Time (s)'), ylabel('Temperature [K]')
set(gcf,'Color','w')

%------------------------------------------------------------------------------------------
% Volume Plot
%------------------------------------------------------------------------------------------
subplot(1,2,2)
plot(t,Vh_temp',t,Vc_temp','-.')
legend('1','2','3','4','5','6','7','8','9','1','2','3','4','5','6','7','8','9')
xlim([0 max(t)]), xlabel('Time (s)'), ylabel('Volume [m^3]')
set(gcf,'Color','w')

%------------------------------------------------------------------------------------------
% All Separate Plots
%------------------------------------------------------------------------------------------
figure
subplot(2,2,1), title('Hot Tank Temperature')
plot(t,Th_temp')
legend('1','2','3','4','5','6','7','8','9')
xlim([0 max(t)]), xlabel('Time (s)'), ylabel('Temperature [K]')
set(gcf,'Color','w')

subplot(2,2,2), title('Hot Tank Volume')
plot(t,Vh_temp')
legend('1','2','3','4','5','6','7','8','9')
xlim([0 max(t)]), xlabel('Time (s)'), ylabel('Volume [m^3]')
set(gcf,'Color','w')

subplot(2,2,3), title('Cold Tank Temperature')
plot(t,Tc_temp','-.')
legend('1','2','3','4','5','6','7','8','9')
xlim([0 max(t)]), xlabel('Time (s)'), ylabel('Temperature [K]')
set(gcf,'Color','w')

subplot(2,2,4), title('Cold Tank Volume')
plot(t,Vc_temp','-.')
legend('1','2','3','4','5','6','7','8','9')
xlim([0 max(t)]), xlabel('Time (s)'), ylabel('Volume [m^3]')
set(gcf,'Color','w')


