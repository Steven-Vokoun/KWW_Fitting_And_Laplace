% Request data from oscilloscope with NIVISA over LAN
[timeOut, dataOut, ~] = acquireOscilloscopeData('10.0.0.3', 1); % Calls acquireOscilloscopeData.m and determineAcquisitionSettings.m

%% Variables
Gain = 1;
Vin = .1;
TestResistor = 931;

%% Setup Data
Current = (dataOut * Gain) / TestResistor;
L = length(Current);
%Normalization_Offset = min(Current);
[MaxCurrent, i0] = max(Current);
Normalization_Offset = 0;
Current = Current - Normalization_Offset;

MinCurrent = mean(Current(1:100)); % Average first 100 samples
Iinf = mean(Current(L-100:L));  % Average last 100 samples

%% Obtain Proper Time Base
TimeCalc = timeOut(i0:end);
CurrentCalc = Current(i0:end);
VoltageCalc = CurrentCalc * TestResistor;

minTime = min(TimeCalc);
TimeCalc = TimeCalc - minTime;

%% Ohm's Law Calculation
RS = Vin / MaxCurrent;
RP = (Vin / (Iinf - MinCurrent)) - RS;
It = (.367 * (MaxCurrent - Iinf)) + Iinf;
[~, ix] = min(abs(CurrentCalc - It));
t = TimeCalc(ix);
Ceff_Randles_Calculated = (t * (RS + RP)) / (RS * RP);
disp(Ceff_Randles_Calculated);
Ohms = ((Vin / RS) * exp(-TimeCalc/t)) + ((Vin / (RS + RP)) * (1 - exp(-TimeCalc/t))) + MinCurrent;
r = corrcoef(Ohms, CurrentCalc);
R2 = r .* r;
Ohms_R2 = R2(1, 2);
disp(RS)
disp(RP)
RP = abs(RP);
% Randles Fitting
try
    [xRData, yRData] = prepareCurveData(TimeCalc, CurrentCalc);
    % Set up fittype and options.
    ft = fittype('((V/S)*exp(-x/((S*P*C)/(S+P)))) + ((V/(S+P))*(1-exp(-x/((S*P*C)/(S+P))))) + M', 'independent', 'x', 'dependent', 'y');
    opts = fitoptions('Method', 'NonlinearLeastSquares');
    opts.Display = 'Off';
    opts.DiffMinChange = 1e-25;
    opts.DiffMaxChange = 1E-15;
    opts.Lower = [1E-10, -.001, RP, RS, Vin];
    opts.MaxFunEvals = 10000;
    opts.MaxIter = 10000;
    opts.StartPoint = [Ceff_Randles_Calculated, MinCurrent, RP, RS, Vin];
    opts.TolFun = 1e-20;
    opts.TolX = 1e-20;
    opts.Upper = [5E-5, .001, RP, RS, Vin];
    
    % Fit model to data.
    [Randlesresult, gof] = fit(xRData, yRData, ft, opts);
    Ceff_Randles_FIT = Randlesresult.C;
    Ceff_Randles_R2 = gof.rsquare;
    disp(Ceff_Randles_FIT);
    disp(Ceff_Randles_R2);
    Randlesdata = feval(Randlesresult, TimeCalc);
catch
    Ceff_Randles_FIT = 0;
    Ceff_Randles_R2 = 0;
    ERROR = 1;
end

%% Randles Fitting with N Value
try
    [xData, yData] = prepareCurveData(TimeCalc, CurrentCalc);
    % Set up fittype and options.
    ft = fittype('((V/S)*exp(-(x/((S*P*C)/(S+P)))^N)) + ((V/(S+P))*(1-(exp(-(x/((S*P*C)/(S+P)))^N)))) + M', 'independent', 'x', 'dependent', 'y');
    opts = fitoptions('Method', 'NonlinearLeastSquares');
    opts.DiffMinChange = 1e-25;
    opts.DiffMaxChange = 1E-15;
    opts.Display = 'Off';
    opts.Lower = [5E-10, -.001, .05, RP, RS, Vin];
    opts.MaxFunEvals = 1000;
    opts.MaxIter = 1000;
    opts.StartPoint = [Ceff_Randles_Calculated, MinCurrent, .4, RP, RS, Vin];
    opts.TolFun = 1e-20;
    opts.TolX = 1e-20;
    opts.Upper = [5E-7, .001, .95, RP, RS, Vin];
    
    % Fit model to data.
    [Nresult, Ngof] = fit(xData, yData, ft, opts);
    Ceff_Randles_FIT_N = Nresult.C;
    KWW_N = Nresult.N;
    Ceff_Randles_R2_N = Ngof.rsquare;
    disp(Ceff_Randles_FIT_N);
    disp(KWW_N);
    disp(Ceff_Randles_R2_N);
    Ndata = feval(Nresult, TimeCalc);
catch
    Ceff_Randles_FIT_N = 0;
    KWW_N = 0;
    Ceff_Randles_R2_N = 0;
    ERROR = 1;
    Ndata = 0;
end




% %% Plotting
% 
% % Create a figure and make it visible
% Combined = figure('visible', 'on');
% 
% % First subplot
% subplot(1,2,1)
% hold on;
% set(gca, 'FontSize', 14);
% 
% % Plot data
% plot(TimeCalc(1:50:end)/1E-3, CurrentCalc(1:50:end)/1E-6, 'Color', 'K', 'LineStyle', 'none', 'Marker', '.', 'LineWidth', 2); % Experimental
% plot(TimeCalc(1:50:end)/1E-3, Randlesdata(1:50:end)/1E-6, 'Color', [0.4 0.4 0.4], 'LineWidth', 2, 'LineStyle', '--'); % Randles
% plot(TimeCalc(1:50:end)/1E-3, Ndata(1:50:end)/1E-6, 'Color', [0.8 0.8 0.8], 'LineWidth', 2); % Randles with N
% 
% % Legend for the first subplot
% lgd1 = legend('Raw Experimental Time Data', ...
%     ['Real Capacitor:  R^2 = ' num2str(Ceff_Randles_R2, 4)], ...
%     ['KWW Function:  R^2 = ' num2str(Ceff_Randles_R2_N, 4)], 'FontSize', 12, 'Location', 'south');
% 
% grid on
% xlabel('Time (ms)', 'FontSize', 14);
% ylabel('Current µA', 'FontSize', 14);
% title("A. Linear Scale", 'FontSize', 16);
% xlim([TimeCalc(1)/1E-3 TimeCalc(end)/1E-3])
% axis square
% 
% % Second subplot
% subplot(1,2,2)
% hold on;
% set(gca, 'XScale', 'log');
% set(gca, 'FontSize', 14);
% 
% % Plot data
% plot(TimeCalc(1:50:end)/1E-3, CurrentCalc(1:50:end)/1E-6, 'Color', 'K', 'LineStyle', 'none', 'Marker', '.', 'LineWidth', 1); % Experimental 
% plot(TimeCalc(1:50:end)/1E-3, Randlesdata(1:50:end)/1E-6, 'Color', [0.4 0.4 0.4], 'LineWidth', 2, 'LineStyle', '--'); % Randles
% plot(TimeCalc(1:50:end)/1E-3, Ndata(1:50:end)/1E-6, 'Color', [0.8 0.8 0.8], 'LineWidth', 2); % Randles with N
% 
% % Legend for the second subplot
% lgd2 = legend('Raw Experimental Time Data', ...
%     ['Real Capacitor:  R^2 = ' num2str(Ceff_Randles_R2, 4)], ...
%     ['KWW Function:  R^2 = ' num2str(Ceff_Randles_R2_N, 4)], 'FontSize', 12, 'Location', 'south');
% 
% grid on
% xlabel('Time (ms)', 'FontSize', 14);
% ylabel('Current µA', 'FontSize', 14);
% title("B. Log Scale", 'FontSize', 14);
% 
% axis square
% 
% % Adjust the legend position for the entire figure
% lgd1.Position(1) = 0.5 - lgd1.Position(3)/2; % Center the legend horizontally
% lgd1.Position(2) = 0.05; % Adjust the vertical position to be at the bottom
% 
% lgd2.Visible = 'off'; % Hide the legend for the second subplot










%% Plot Graph
% Combined = figure('visible', 'on');
% subplot(1,2,1)
% hold on;
% set(gca, 'FontSize', 14);
% plot(TimeCalc(1:50:end), CurrentCalc(1:50:end)/1E-6, 'Color', 'K', 'LineStyle', 'none', 'Marker', '.', 'LineWidth', 2); % Experimental
% plot(TimeCalc(1:50:end), Randlesdata(1:50:end)/1E-6, 'Color', [86/255, 180/255, 233/255], 'LineWidth', 2); % Randles
% plot(TimeCalc(1:50:end), Ndata(1:50:end)/1E-6, 'Color', [0/255, 158/255, 115/255], 'LineStyle', '--', 'LineWidth', 2); % Randles with N
% 
% legend('Raw Experimental Time Data', ...
%     ['Real Capacitor:  R^2 = ' num2str(Ceff_Randles_R2, 4)], ...
%     ['KWW Function:  R^2 = ' num2str(Ceff_Randles_R2_N, 4)], 'FontSize', 12)
% 
% grid on
% xlabel('Time (s)', 'FontSize', 14);
% ylabel('Current µA', 'FontSize', 14);
% title("A. Linear Scale", 'FontSize', 16);
% 
% axis square
% 
% subplot(1,2,2)
% hold on;
% set(gca, 'XScale', 'log');
% set(gca, 'FontSize', 14);
% plot(TimeCalc(1:50:end), CurrentCalc(1:50:end)/1E-6, 'Color', 'K', 'LineStyle', 'none', 'Marker', '.', 'LineWidth', 2); % Experimental 
% plot(TimeCalc(1:50:end), Randlesdata(1:50:end)/1E-6, 'Color', [86/255, 180/255, 233/255], 'LineWidth', 2); % Randles
% plot(TimeCalc(1:50:end), Ndata(1:50:end)/1E-6, 'Color', [0/255, 158/255, 115/255], 'LineStyle', '--', 'LineWidth', 2); % Randles with N
% 
% legend('Raw Experimental Time Data', ...
%     ['Real Capacitor:  R^2 = ' num2str(Ceff_Randles_R2, 4)], ...
%     ['KWW Function:  R^2 = ' num2str(Ceff_Randles_R2_N, 4)], 'FontSize', 12)
% 
% grid on
% xlabel('Time (s)', 'FontSize', 14);
% ylabel('Current µA', 'FontSize', 14);
% title("B. Log Scale", 'FontSize', 14);
% 
% axis square
% 
% RS = RS - TestResistor;
