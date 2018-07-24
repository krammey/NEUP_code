% addpath('/Users/trins/Documents/MATLAB/altmany-export_fig-9d97e2c')
% addpath('/Users/trins/Documents/MATLAB/mtit')
% 
% start_time = tic;
% Revenue_Final = [];
% annualCC_Final = [];
% TotalCC_Final = [];
% annualOM_Final = [];
% fixedOM_Final = [];
% for caseNumber = 0:6
%     [Revenue_Results,annualCC_Results,annualOM_Results,TotalCC_Results,fixedOM_Results] = run_steam_model(caseNumber);
%     nrows          = size(Revenue_Results,1);
%     Revenue_Final  = [Revenue_Final; caseNumber*ones(nrows,1) Revenue_Results];
%     annualCC_Final = [annualCC_Final; caseNumber*ones(nrows,1) annualCC_Results];
%     TotalCC_Final  = [TotalCC_Final; caseNumber*ones(nrows,1) TotalCC_Results];
%     annualOM_Final = [annualOM_Final; caseNumber*ones(nrows,1) annualOM_Results];
%     fixedOM_Final  = [fixedOM_Final; caseNumber*ones(nrows,1) fixedOM_Results];
% end
% fprintf('Total run time = %.2f seconds.\n', toc(start_time));
% disp(' ')
% 
% %============================================
% % Plotting
% %============================================
% T = size(Revenue_Final,2)-2; % length of T_MAX_array in run_steam_model.m
% red = "0.85 0.33 0.1"; blue = "0 0.45 0.74"; grey = "0.4 0.4 0.4";
% MARKERS = ["-d";"-s";"-o"];
% COLORS = [blue; red; grey];
% LABELS = ["Annual Net Revenue [$/kWh_e]"; "Annual CC [$/kWh_e]"; "Total O&M [$/kWh_e]"; "Total CC [$/kWh_e]";"Fixed O&M [$/kW-year]"]; % new labels after converting units
% WIDTH = 1;
% BaseCase = 0;
% 
% % Re-organize results to be able to plot using for loop 
% Results        = [];
% Results(:,:,1) = Revenue_Final;
% Results(:,:,2) = annualCC_Final;
% Results(:,:,3) = annualOM_Final;
% Results(:,:,4) = TotalCC_Final;
% Results(:,:,5) = fixedOM_Final;
% 
% % Extract power array
% POWER_array    = unique(Results(:,2,:));
% P              = length(POWER_array);

figure
for r=1:5 % loop through results: 1-Revenue, 2-CC, 3-Total OM, 4-Total CC, 5-Fixed OM Results
    subplot(2,3,r)
    data_temp = Results(:,:,r);
    hold on
    for idx = 1:P
        PlotData_temp = data_temp(data_temp(:,2)==POWER_array(idx),:);
        plot(categorical(1:T),PlotData_temp(:,3:T+2),':','Color',COLORS(idx),'LineWidth',WIDTH)
        plot(categorical(1:T),PlotData_temp(4,3:T+2),'-','Color',COLORS(idx),'LineWidth',1)
    end
    ylabel(LABELS(r)), xlabel('Hours of storage')
    hold off
end

subplot(2,3,5)
hold on
p1 = plot(0,'-','Color',COLORS(1));
p2 = plot(0,'-','Color',COLORS(2));
p3 = plot(0,'-','Color',COLORS(3));
p4 = plot(0,':','Color',COLORS(2));
p5 = plot(0,':','Color',COLORS(1));
p6 = plot(0,':','Color',COLORS(3));
hold off
LEGEND = ["50 MW [Base Case]" "100 MW [Base Case]" "150 MW [Base Case]" "50 MW [All Other Cases]" "100 MW [All Other Cases]" "150 MW [All Other Cases]"];
legend([p1 p2 p3 p4 p5 p6],LEGEND,'Location','BestOutside'),legend('boxoff')
set(gcf,'Color','w')
set(findall(gcf,'-property','FontSize'),'FontSize',10)
set(findall(gcf,'-property','MarkerSize'),'MarkerSize',4)

% export_fig(['Case' num2str(caseNumber) '.fig'])
disp(' ')
% 
% 
% 
% 
% 
% %----------------------
% % Net Revenue
% %----------------------
% subplot(1,5,1)
% hold on
% % 50 MW - Blue
% Revenue_50 = Revenue_Final(Revenue_Final(:,2)==50,:);
% plot(categorical(1:T),Revenue_50(:,3:10),':','Color',blue,'LineWidth',.5)
% plot(categorical(1:T),Revenue_50(4,3:10),'- .','MarkerFaceColor',blue,'Color',blue,'LineWidth',1)
% % 100 MW - Red
% Revenue_100 = Revenue_Final(Revenue_Final(:,2)==100,:);
% plot(categorical(1:T),Revenue_100(:,3:10),':','Color',red,'LineWidth',.5)
% plot(categorical(1:T),Revenue_100(4,3:10),'- .','MarkerFaceColor',red,'Color',red,'LineWidth',1)
% % 150 MW - Yellow
% Revenue_150 = Revenue_Final(Revenue_Final(:,2)==150,:);
% plot(categorical(1:T),Revenue_150(:,3:10),':','Color',grey,'LineWidth',.5)
% plot(categorical(1:T),Revenue_150(4,3:10),'- .','MarkerFaceColor',grey,'Color',grey,'LineWidth',1)
% ylabel("Annual Net Revenue [million $]"), xlabel('Hours of storage')
% hold off
% 
% %----------------------
% % CC
% %----------------------
% subplot(1,5,2)
% hold on
% % 50 MW - Blue
% CC_50 = annualCC_Final(annualCC_Final(:,2)==50,:);
% plot(categorical(1:T),CC_50(:,3:10),':','Color',blue,'LineWidth',WIDTH)
% plot(categorical(1:T),CC_50(4,3:10),'- .','MarkerFaceColor',blue,'Color',blue,'LineWidth',1)
% % 100 MW - Red
% CC_100 = annualCC_Final(annualCC_Final(:,2)==100,:);
% plot(categorical(1:T),CC_100(:,3:10),':','Color',red,'LineWidth',WIDTH)
% plot(categorical(1:T),CC_100(4,3:10),'- .','MarkerFaceColor',red,'Color',red,'LineWidth',1)
% % 150 MW - Yellow
% CC_150 = annualCC_Final(annualCC_Final(:,2)==150,:);
% plot(categorical(1:T),CC_150(:,3:10),':','Color',grey,'LineWidth',WIDTH)
% plot(categorical(1:T),CC_150(4,3:10),'- .','MarkerFaceColor',grey,'Color',grey,'LineWidth',1)
% ylabel("Annual CC [million $]"), xlabel('Hours of storage')
% hold off
% %----------------------
% % Total OM
% %----------------------
% subplot(1,5,3)
% hold on
% % 50 MW - Blue
% TotalOM_50 = annualOM_Final(annualOM_Final(:,2)==50,:);
% plot(categorical(1:T),TotalOM_50(:,3:10),':','Color',blue,'LineWidth',WIDTH)
% plot(categorical(1:T),TotalOM_50(4,3:10),'- .','MarkerFaceColor',blue,'Color',blue,'LineWidth',1)
% % 100 MW - Red
% TotalOM_100 = annualOM_Final(annualOM_Final(:,2)==100,:);
% plot(categorical(1:T),TotalOM_100(:,3:10),':','Color',red,'LineWidth',WIDTH)
% plot(categorical(1:T),TotalOM_100(4,3:10),'- .','MarkerFaceColor',red,'Color',red,'LineWidth',1)
% % 150 MW - Yellow
% TotalOM_150 = annualOM_Final(annualOM_Final(:,2)==150,:);
% plot(categorical(1:T),TotalOM_150(:,3:10),':','Color',grey,'LineWidth',WIDTH)
% plot(categorical(1:T),TotalOM_150(4,3:10),'- .','MarkerFaceColor',grey,'Color',grey,'LineWidth',1)
% ylabel('Total O&M [million $]'), xlabel('Hours of storage')
% hold off
% %----------------------
% % Total CC
% %----------------------
% subplot(1,5,4)
% hold on
% % 50 MW - Blue
% TotalCC_50 = TotalCC_Final(TotalCC_Final(:,2)==50,:);
% Pdot1 = plot(categorical(1:T),TotalCC_50(:,3:10),':','Color',blue,'LineWidth',WIDTH);
% p1 = plot(categorical(1:T),TotalCC_50(4,3:10),'- .','MarkerFaceColor',blue,'Color',blue,'LineWidth',1);
% % 100 MW - Red
% TotalCC_100 = TotalCC_Final(TotalCC_Final(:,2)==100,:);
% Pdot2 = plot(categorical(1:T),TotalCC_100(:,3:10),':','Color',red,'LineWidth',WIDTH);
% p2 = plot(categorical(1:T),TotalCC_100(4,3:10),'- .','MarkerFaceColor',red,'Color',red,'LineWidth',1);
% % 150 MW - Yellow
% TotalCC_150 = TotalCC_Final(TotalCC_Final(:,2)==150,:);
% Pdot3 = plot(categorical(1:T),TotalCC_150(:,3:10),':','Color',grey,'LineWidth',WIDTH);
% p3 = plot(categorical(1:T),TotalCC_150(4,3:10),'- .','MarkerFaceColor',grey,'Color',grey,'LineWidth',1);
% 
% 
% legend([p1 p2 p3 Pdot1(1) Pdot2(1) Pdot3(1)],["50 MW [Base Case]" "100 MW [Base Case]" "150 MW [Base Case]" "50 MW [All Other Cases]" "100 MW [All Other Cases]" "150 MW [All Other Cases]"],'Location','BestOutside'),legend('boxoff')
% hold off
% ylabel('Total CC [million $]'), xlabel('Hours of storage')
% %----------------------
% set(gcf,'Color','w')
% set(findall(gcf,'-property','FontSize'),'FontSize',10)
% set(findall(gcf,'-property','MarkerSize'),'MarkerSize',4)
% %----------------------
% 
% 
% 
