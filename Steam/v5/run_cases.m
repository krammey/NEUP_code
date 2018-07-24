addpath('/Users/trins/Documents/MATLAB/altmany-export_fig-9d97e2c')
addpath('/Users/trins/Documents/MATLAB/mtit')

start_time = tic;
Revenue_Final = [];
CC_Final = [];
TotalCC_Final = [];
TotalOM_Final = [];
for caseNumber = 0:6
    [Revenue_Results,CC_Results,TotalOM_Results,TotalCC_Results] = run_steam_model(caseNumber);
    nrows = size(Revenue_Results,1);
    Revenue_Final = [Revenue_Final; caseNumber*ones(nrows,1) Revenue_Results];
    CC_Final = [CC_Final; caseNumber*ones(nrows,1) CC_Results];
    TotalCC_Final = [TotalCC_Final; caseNumber*ones(nrows,1) TotalCC_Results];
    TotalOM_Final = [TotalOM_Final; caseNumber*ones(nrows,1) TotalOM_Results];
end
fprintf('Total run time = %.2f seconds.\n', toc(start_time));
disp(' ')

%----------------------
% Save Results
%----------------------
Final_Results = [Revenue_Final; zeros(1,size(Revenue_Final,2));...
                 CC_Final;      zeros(1,size(CC_Final,2));...
                 TotalCC_Final; zeros(1,size(TotalCC_Final,2));...
                 TotalOM_Final; zeros(1,size(TotalOM_Final,2))];
% csvwrite('All_Cases_Results.csv', Final_Results)

%======================
% Plotting
%======================
T = size(Revenue_Final,2)-2;
red = '0.85 0.33 0.1';
blue = '0 0.45 0.74';
grey = '0.4 0.4 0.4';
% yellow = '0.93 0.69 0.13';
% purple = '0.6350, 0.0780, 0.1840';
WIDTH = 1;
figure
%----------------------
% Net Revenue
%----------------------
subplot(1,5,1)
hold on
% 50 MW - Blue
Revenue_50 = Revenue_Final(Revenue_Final(:,2)==50,:);
plot(categorical(1:T),Revenue_50(:,3:10),':','Color',blue,'LineWidth',.5)
plot(categorical(1:T),Revenue_50(4,3:10),'- .','MarkerFaceColor',blue,'Color',blue,'LineWidth',1)
% 100 MW - Red
Revenue_100 = Revenue_Final(Revenue_Final(:,2)==100,:);
plot(categorical(1:T),Revenue_100(:,3:10),':','Color',red,'LineWidth',.5)
plot(categorical(1:T),Revenue_100(4,3:10),'- .','MarkerFaceColor',red,'Color',red,'LineWidth',1)
% 150 MW - Yellow
Revenue_150 = Revenue_Final(Revenue_Final(:,2)==150,:);
plot(categorical(1:T),Revenue_150(:,3:10),':','Color',grey,'LineWidth',.5)
plot(categorical(1:T),Revenue_150(4,3:10),'- .','MarkerFaceColor',grey,'Color',grey,'LineWidth',1)
ylabel("Annual Net Revenue [million $]"), xlabel('Hours of storage')
hold off
%----------------------
% CC
%----------------------
subplot(1,5,2)
hold on
% 50 MW - Blue
CC_50 = CC_Final(CC_Final(:,2)==50,:);
plot(categorical(1:T),CC_50(:,3:10),':','Color',blue,'LineWidth',WIDTH)
plot(categorical(1:T),CC_50(4,3:10),'- .','MarkerFaceColor',blue,'Color',blue,'LineWidth',1)
% 100 MW - Red
CC_100 = CC_Final(CC_Final(:,2)==100,:);
plot(categorical(1:T),CC_100(:,3:10),':','Color',red,'LineWidth',WIDTH)
plot(categorical(1:T),CC_100(4,3:10),'- .','MarkerFaceColor',red,'Color',red,'LineWidth',1)
% 150 MW - Yellow
CC_150 = CC_Final(CC_Final(:,2)==150,:);
plot(categorical(1:T),CC_150(:,3:10),':','Color',grey,'LineWidth',WIDTH)
plot(categorical(1:T),CC_150(4,3:10),'- .','MarkerFaceColor',grey,'Color',grey,'LineWidth',1)
ylabel("Annual CC [million $]"), xlabel('Hours of storage')
hold off
%----------------------
% Total OM
%----------------------
subplot(1,5,3)
hold on
% 50 MW - Blue
TotalOM_50 = TotalOM_Final(TotalOM_Final(:,2)==50,:);
plot(categorical(1:T),TotalOM_50(:,3:10),':','Color',blue,'LineWidth',WIDTH)
plot(categorical(1:T),TotalOM_50(4,3:10),'- .','MarkerFaceColor',blue,'Color',blue,'LineWidth',1)
% 100 MW - Red
TotalOM_100 = TotalOM_Final(TotalOM_Final(:,2)==100,:);
plot(categorical(1:T),TotalOM_100(:,3:10),':','Color',red,'LineWidth',WIDTH)
plot(categorical(1:T),TotalOM_100(4,3:10),'- .','MarkerFaceColor',red,'Color',red,'LineWidth',1)
% 150 MW - Yellow
TotalOM_150 = TotalOM_Final(TotalOM_Final(:,2)==150,:);
plot(categorical(1:T),TotalOM_150(:,3:10),':','Color',grey,'LineWidth',WIDTH)
plot(categorical(1:T),TotalOM_150(4,3:10),'- .','MarkerFaceColor',grey,'Color',grey,'LineWidth',1)
ylabel('Total O&M [million $]'), xlabel('Hours of storage')
hold off
%----------------------
% Total CC
%----------------------
subplot(1,5,4)
hold on
% 50 MW - Blue
TotalCC_50 = TotalCC_Final(TotalCC_Final(:,2)==50,:);
Pdot1 = plot(categorical(1:T),TotalCC_50(:,3:10),':','Color',blue,'LineWidth',WIDTH);
p1 = plot(categorical(1:T),TotalCC_50(4,3:10),'- .','MarkerFaceColor',blue,'Color',blue,'LineWidth',1);
% 100 MW - Red
TotalCC_100 = TotalCC_Final(TotalCC_Final(:,2)==100,:);
Pdot2 = plot(categorical(1:T),TotalCC_100(:,3:10),':','Color',red,'LineWidth',WIDTH);
p2 = plot(categorical(1:T),TotalCC_100(4,3:10),'- .','MarkerFaceColor',red,'Color',red,'LineWidth',1);
% 150 MW - Yellow
TotalCC_150 = TotalCC_Final(TotalCC_Final(:,2)==150,:);
Pdot3 = plot(categorical(1:T),TotalCC_150(:,3:10),':','Color',grey,'LineWidth',WIDTH);
p3 = plot(categorical(1:T),TotalCC_150(4,3:10),'- .','MarkerFaceColor',grey,'Color',grey,'LineWidth',1);
p4 = plot(0,'-k');
p5 = plot(0,':k');
legend([p1 p2 p3 Pdot1(1) Pdot2(1) Pdot3(1)],["50 MW [Base Case]" "100 MW [Base Case]" "150 MW [Base Case]" "50 MW [All Other Cases]" "100 MW [All Other Cases]" "150 MW [All Other Cases]"],'Location','BestOutside'),legend('boxoff')
hold off
ylabel('Total CC [million $]'), xlabel('Hours of storage')
%----------------------
set(gcf,'Color','w')
set(findall(gcf,'-property','FontSize'),'FontSize',10)
set(findall(gcf,'-property','MarkerSize'),'MarkerSize',4)
%----------------------



