function [Th,Tc,Vh,Vc] = run_salt_tank(HXC,T_HL,T_CL,SALT)
% This function out puts 4 arrays that depict the parameter value over
% time. The arrays are for hot and cold tanks temps and hot and cold tank
% volumes, respoectively. Inputs are only those that will be varied in the
% sensitivity analysis.

%==========================================================================================
% Inputs that are not changed in sensitivity analysis
%==========================================================================================
% Salt Tank Inputs  | Value             | Units       | Description/Notes
%-------------------|-------------------|-------------|-----------------------------
radius              = 6;                % m             both tanks have same dimensions
height              = 6;                % m             in Milestone 3 report a tank height of 12 was used
init_hot_vol_frac   = 0.75;             % -             total fluid volume equals the volume of one tank.  This is the fraction of the total fluid that's initially in the hot tank
P_in                = 50E6;             % J/s           rate at which THERMAL power is transferred in when charging
P_out               = 30E6;             % J/s           rate at which THERMAL power is transferred out when discharging
T_env               = 300;              % K             ambient temperature outside the tank
if SALT == "Solar"
    C_p             = 1496;             % J/kg K        heat capacity of salt
    rho             = 1899;             % kg/m3         density of salt
elseif SALT == "Hitec XL"
    C_p             = 1447;
    rho             = 1992;
elseif SALT == "Glauber"
    C_p             = 3349;
    rho             = 1460;
else
    print("Error in salt input");
end

%-------------------------------------------------------------------------
% Initial Calculations
%-------------------------------------------------------------------------
volume                  = pi*radius*radius*height;       % m3, tank volume
side_area               = 2*pi*radius*height;            % m2, area of tank sides


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
event(4)                    = 1;                % do nothing for a while
event(5)                    = 0;                % discharge
event(6)                    = 2;                % charge
event(7)                    = 1;                % do nothing for a while
event(8)                    = 0;                % discharge    
endt_event(1)               = 1800;             % s             event 1 ends after 1800 s (30 minutes)
endt_event(2)               = 3600;             % s             event 2 ends at 3600 s (hence lasts 1800 s)
endt_event(3)               = 5400;             % s             sets the end time of each event
endt_event(4)               = 12000;            % s
endt_event(5)               = 15000;            % s
endt_event(6)               = 17600;            % s
endt_event(7)               = 36000;            % s
endt_event(8)               = 38400;            % s
dt                          = 100;              % s             time step in seconds - I haven't tested convergence, could possibly jack this way up
tsteps_per_event(1)     = round(endt_event(1)/dt);       % number of time steps in the first event
endt_event(1)           = dt*tsteps_per_event(1);
nt                      = tsteps_per_event(1);
for i=2:N_EVENTS        % calculates total number of time steps for all events
    tsteps_per_event(i) = round((endt_event(i)-endt_event(i-1))/dt);% number of time steps in event i
    endt_event(i)       = endt_event(i-1) + dt*tsteps_per_event(i); % rounds entries here?
    nt                  = nt + tsteps_per_event(i);
end
t                       = (0:dt:(nt-1)*dt); % time array?

%==========================================================================================
% end of inputs
%==========================================================================================


%-------------------------------------------------------------------------
% Creates arrays of length nt for time, temperatures, volumes, and masses
%-------------------------------------------------------------------------
Th                      = [1,nt];
Tc                      = [1,nt];
Vh                      = [1,nt];
Vc                      = [1,nt];
mh                      = [1,nt];
mc                      = [1,nt];
% mdotchg = zeros(nt,1)';
% mdotdis = zeros(nt,1)';


Th(1) = T_HL;                         % initial T in hot tank
Tc(1) = T_CL;                         % initial T in cold tank
Vh(1) = volume*init_hot_vol_frac;     % initial volume in hot tank
Vc(1) = volume*(1-init_hot_vol_frac); % initial volume in cold tank
% mdotchg(1) = 0;
% mdotdis(1) = 0;
ix = 2;
% mt = [];

% calculated temperatures and volumes at the time steps for each case
% the formulas used for Th, Tc, Vh, and Vc can be found in the report -- derived by Schneider 
for i=1:N_EVENTS
    if      event(i) == 0  % discharge
        for j=1:tsteps_per_event(i) 
            Th(ix) = Th(ix-1)-dt*HXC*side_area*(Th(ix-1)-T_env)/(C_p*rho*Vh(ix-1));
            Tc(ix) = Tc(ix-1)-dt*HXC*side_area*(Tc(ix-1)-T_env)/(C_p*rho*Vc(ix-1));
            Vh(ix) = Vh(ix-1)-dt*P_out/(rho*C_p*(Th(ix-1)-T_CL)); 
            Vc(ix) = Vc(ix-1)+dt*P_out/(rho*C_p*(Th(ix-1)-T_CL));
            % mt     = dt*(P_out/(rho*C_p*(Th(ix-1)-T_CL))); % Why is this the same as the second term in V eqn? Also this doesn't appear to be used later********
            ix     = ix+1;
        end % for int j
        ix = ix-2; % reset time step index after each event
    elseif   event(i) == 1  % just heat loss
        for j=1:tsteps_per_event(i)
            Th(ix) = Th(ix-1)-dt*HXC*side_area*(Th(ix-1)-T_env)/(C_p*rho*Vh(ix-1));
            Tc(ix) = Tc(ix-1)-dt*HXC*side_area*(Tc(ix-1)-T_env)/(C_p*rho*Vc(ix-1));
            Vh(ix) = Vh(ix-1);
            Vc(ix) = Vc(ix-1);
            ix     = ix+1;
        end %for int j
        ix = ix-2; % reset time step index after each event
    elseif  event(i) == 2  % charge
        for j=1:tsteps_per_event(i)
            Th(ix) = Th(ix-1) + dt*( P_in*(T_HL-Th(ix-1))/(C_p*rho*Vh(ix-1)*(T_HL-Tc(ix-1))) - HXC*side_area*(Th(ix-1)-T_env)/(C_p*rho*Vh(ix-1)));
            Tc(ix) = Tc(ix-1) - dt*( P_in*(T_HL-Th(ix-1))/(C_p*rho*Vh(ix-1)*(T_HL-Tc(ix-1))) + HXC*side_area*(Tc(ix-1)-T_env)/(C_p*rho*Vc(ix-1)));
            Vh(ix) = Vh(ix-1) + dt*(P_in/(rho*C_p*(T_HL-Tc(ix-1))));
            Vc(ix) = Vc(ix-1) - dt*(P_in/(rho*C_p*(T_HL-Tc(ix-1)))); 
            ix     = ix+1;
        end %for int j
        ix = ix-2; % reset time step index after each event
    else
        break;
    end
end % for int i




end

