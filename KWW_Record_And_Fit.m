%function [Ceff_Randles_Calculated, Ohms_R2, Ceff_Randles_FIT, Ceff_Randles_R2, Ceff_Randles_FIT_N, Ceff_Randles_R2_N, KWW_N] = KWW_Record_And_Fit(Vin, Title)
%KWW_Record_And_Fit communicates with the python script to record ohms law,
%randles, and KWW fitting of platinum in PBS
%   Dependancy on acquireOscilloscopeData.m and
%   determineAcquisitionSettings.m
%   https://github.com/dadul96/Siglent_SDS1202X-E_VISA_MATLAB
%% Request data from oscilloscope with NIVISA over LAN
[timeOut, dataOut, ~] = acquireOscilloscopeData('192.168.1.3',1);  %% Calls accquireOscilloscopeData.m and determineAcquisitionSettings.m

%% Variables
%Vin = .1;    %takes in from the function call
%Title = -.9;
Offset = .1;  % time out
TestResistor = 100000; %100k

%% Setup Data
timeOut = timeOut + Offset;
Current = dataOut / TestResistor;
[MaxCurrent,i0] = max(Current);
L = length (Current);
Normalization_Offset = min(Current);
Current = Current - Normalization_Offset;
MinCurrent = mean(Current(1:100)); %average first 100 samples
Iinf = mean(Current(L-100:L));  %average last 100 samples

%% Obtain Proper Time Base
TimeCalc = timeOut(i0:end);
CurrentCalc = Current(i0:end);
CurrentCalc = movmean(CurrentCalc,50);

%% Ohms Law Calculation
RS = Vin/MaxCurrent;
RP = (Vin/(Iinf-MinCurrent))-RS;
It = (.367*(MaxCurrent-Iinf))+Iinf;
[ ~ , ix ] = min( abs(CurrentCalc-It));
t = TimeCalc(ix);
Ceff_Randles_Calculated = (t*(RS+RP))/(RS*RP);
disp(Ceff_Randles_Calculated);
Ohms = ((Vin/RS)*exp(-TimeCalc/t))+((Vin/(RS+RP))*(1-exp(-TimeCalc/t)))+MinCurrent;
r = corrcoef(Ohms, CurrentCalc);
R2 = r.*r;
Ohms_R2 = R2(1,2);

%% Randles Fitting
[xRData, yRData] = prepareCurveData(TimeCalc, CurrentCalc);
% Set up fittype and options.
ft = fittype('((V/S)*exp(-x/(S*C))) + ((V/(S+P))*(1-exp(-x/(S*C)))) + M', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions('Method','NonlinearLeastSquares');
opts.Display = 'Off';
opts.DiffMinChange = 1e-25;
opts.DiffMaxChange = 1E-15;
opts.Lower = [1E-12 -.001 RP RS Vin];
opts.MaxFunEvals = 100000;
opts.MaxIter = 100000;
opts.StartPoint = [Ceff_Randles_Calculated MinCurrent RP RS Vin];
opts.TolFun = 1e-20;
opts.TolX = 1e-20;
opts.Upper = [1E-6 .001 RP RS Vin];

% Fit model to data.
[Randlesresult, gof] = fit( xRData, yRData, ft, opts );
Ceff_Randles_FIT = Randlesresult.C;
Ceff_Randles_R2 = gof.rsquare;
disp(Ceff_Randles_FIT);
disp(Ceff_Randles_R2);
Randlesdata = feval(Randlesresult,TimeCalc);

%% Randles Fitting with N Value
[xData, yData] = prepareCurveData( TimeCalc, CurrentCalc );
% Set up fittype and options.
ft = fittype( '((V/S)*exp(-(x/(S*C))^N)) + ((V/(S+P))*(1-(exp(-(x/(S*C))^N)))) + M', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions('Method','NonlinearLeastSquares');
opts.DiffMinChange = 1e-25;
opts.DiffMaxChange = 1E-15;
opts.Display = 'Off';
opts.Lower = [1E-12 MinCurrent-.001 .01 RP RS Vin];
opts.MaxFunEvals = 100000;
opts.MaxIter = 100000;
opts.StartPoint = [3E-9 MinCurrent .4 RP RS Vin];
opts.TolFun = 1e-20;
opts.TolX = 1e-20;
opts.Upper = [1E-6 MinCurrent+.001 .99 RP RS Vin];


% Fit model to data.
[Nresult, Ngof] = fit( xData, yData, ft, opts );
Ceff_Randles_FIT_N = Nresult.C;
KWW_N = Nresult.N;
Ceff_Randles_R2_N = Ngof.rsquare;
disp(Ceff_Randles_FIT_N);
disp(Ceff_Randles_R2_N);
Ndata = feval(Nresult,TimeCalc);

%% Plot Graph
Combined = figure();
hold on;
%set(gca, 'XScale', 'log');
E = plot(TimeCalc, CurrentCalc, 'Color', 'K', 'LineStyle','-.'); %experimental
O = plot(TimeCalc, Ohms, 'Color', [230/255,159/255,0/255]); %Ohms 
R = plot(TimeCalc, Randlesdata, 'Color', [86/255,180/255,233/255]); %Randles
N = plot(TimeCalc, Ndata, 'Color', [0/255,158/255,115/255], 'LineStyle','--'); %Randles with N



legend('Raw Experimental', ['Ohms Law R^2 = ' num2str(R2(1,2), 4)], ...
    ['Randles R^2 = ' num2str(Ceff_Randles_R2, 4)], ...
    ['KWW Function R^2 = ' num2str(Ceff_Randles_R2_N, 4)])
grid on
xlabel('Time (S)');
ylabel('Current Across Test Resistor');
title(["Transient Response of Platinum Electrode"]);
subtitle([num2str(Title) 'V with a ' num2str(Vin) 'V Step'])

%% Save Raw Data and plots
T = table(TimeCalc, CurrentCalc, Ohms, Randlesdata, Ndata);
writetable(T, [num2str(Title) 'V.txt'], 'Delimiter', '\t');

exportgraphics(Combined, [num2str(Title) 'V.png'], 'Resolution', '300');

%end