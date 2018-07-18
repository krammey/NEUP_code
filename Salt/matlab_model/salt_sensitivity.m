

T_range    = [350 550;...
              375 525;...
              400 500];          % [K] Sets of hot and cold temp ranges. Note temp diff is [200, 150, 100], centered around 450 K
HXC_range  = [20, 35, 50];       % heat transfer coeff for natural circulation convective heat xfer to environment - W/m2/K
SALT_range = ["Solar","Glauber", "Hitec XL"];

[Th,Tc,Vh,Vc] = run_salt_tank(HXC_range(1),T_range(1,2),T_range(1,1),SALT_range(1));

count = 1;
LEGEND = ["Count","HXC","T_HL","T_CL","SALT"];
for k=1:length(SALT_range)
    for i = 1:length(HXC_range)       % i iterates through the heat transfer coefficients
        for j=1:length(T_range)       % j iterates though the temperature ranges
            [Th(count,:),Tc(count,:),Vh(count,:),Vc(count,:)] = run_salt_tank(HXC_range(i),T_range(j,2),T_range(j,1),SALT_range(k));
            LEGEND(count+1,:) = [num2str(count),num2str(HXC_range(i)),num2str(T_range(j,2)),num2str(T_range(j,1)),SALT_range(k)];
            count = count+1;
        end
    end % for int g
end %for int h

time = 100*(1:length(Th));
%==========================================================================================
% Plotting
%==========================================================================================

%------------------------------------------------------------------------------------------
% Temperature Plot
%------------------------------------------------------------------------------------------
figure
subplot(1,2,1)
plot(time,Th',time,Tc','-.'),
xlabel('Time (s)'), ylabel('Temperature [K]')
set(gcf,'Color','w')

%------------------------------------------------------------------------------------------
% Volume Plot
%------------------------------------------------------------------------------------------
subplot(1,2,2)
plot(time,Vh',time,Vc','-.'),
xlabel('Time (s)'), ylabel('Volume [m^3]')
set(gcf,'Color','w')

%------------------------------------------------------------------------------------------
% All Separate Plots
%------------------------------------------------------------------------------------------
figure
subplot(2,2,1), title('Hot Tank Temperature')
plot(Th'), xlabel('Time (s)'), ylabel('Temperature [K]')
set(gcf,'Color','w')

subplot(2,2,2), title('Hot Tank Volume')
plot(Vh'), xlabel('Time (s)'), ylabel('Volume [m^3]')
set(gcf,'Color','w')

subplot(2,2,3), title('Cold Tank Temperature')
plot(Tc','-.'), xlabel('Time (s)'), ylabel('Temperature [K]')
set(gcf,'Color','w')

subplot(2,2,4), title('Cold Tank Volume')
plot(Vc','-.'), xlabel('Time (s)'), ylabel('Volume [m^3]')
set(gcf,'Color','w')


% addpath('/Users/trins/Documents/MATLAB/altmany-export_fig-9d97e2c')
% export_fig('')

