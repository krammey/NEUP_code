
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