clear all

%% Accuire Oscilloscope
[timeOut, dataOut, ~] = acquireOscilloscopeData('10.0.0.3', 1); % Calls acquireOscilloscopeData.m and determineAcquisitionSettings.m

Gain = 1;
Vin = .1;
TestResistor = 9931;

Current = (dataOut * Gain) / TestResistor;
L = length(Current);
Normalization_Offset = min(Current);
Current = Current - Normalization_Offset;
[MaxCurrent, i0] = max(Current);
MinCurrent = mean(Current(1:100)); % Average first 100 samples
Iinf = mean(Current(L-100:L));  % Average last 100 samples

TimeCalc = timeOut(i0:end);
CurrentCalc = Current(i0:end);
VoltageCalc = CurrentCalc * TestResistor;

minTime = min(TimeCalc);
TimeCalc = TimeCalc - minTime;



%% Time to Frequency
V = .1;
t = 0:1E-8:.1;
s = logspace(1,5,100);
GPU_S = gpuArray(s);

Electrode1 = importfile('0(3).txt');


%Real
C = 1.63E-08;
Rp = 2.36E+05;
Rs = 2730;

%KWW
N = .413;
S = 1560;
P = 2.99E+05;
KWW_C = 1.70E-08;

%CPE
a = .741;
Q = 1.11E-07;
RS = 2070;
RP = 9.66E+05;

Time_Function = @(t)V./((((V/Rs)*exp(-t/((Rs*Rp*C)/(Rs+Rp)))) + ((V/(Rs+Rp))*(1-exp(-t/((Rs*Rp*C)/(Rs+Rp)))))));
KWW_Time_Function = @(t)V./(((V./S).*exp(-1.*(t./((S.*P.*KWW_C)./(S+P))).^N)) + ((V./(S+P))*(1-(exp(-(t./((S.*P.*KWW_C)./(S+P))).^N)))));
Time = Time_Function(t);
KWW_Time = KWW_Time_Function(t);

Z_CPE = sqrt((Electrode1.Z1).^2 + (Electrode1.Z2).^2);






LapTime = @(s) LapTrans(Time_Function,s);
LapKWW = @(s) LapTrans(KWW_Time_Function,s);

yTime = arrayfun(LapTime, s);
yKWW = arrayfun(LapKWW, s);


figure()
hold on
set(gca, 'XScale', 'log')
plot(s, yTime)
plot(s, yKWW)
%plot((Electrode1.Freq*(2*pi)), Z_CPE)
title('Time To Frequency Domain')

legend('EIS fit with a Real Capacitor', 'Transient KWW Data', 'EIS as a CPE')





function Data = importfile(filename, dataLines)
%IMPORTFILE Import data from the EISSA File
if nargin < 2
    dataLines = [2, Inf];
end

%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 3);

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = " ";

% Specify column names and types
opts.VariableNames = ["Z1", "Z2", "Freq"];
opts.VariableTypes = ["double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts.ConsecutiveDelimitersRule = "join";
opts.LeadingDelimitersRule = "ignore";

% Import the data
Data = readtable(filename, opts);
end
















% Time = Time_Function(t);
% Time = gpuArray(Time);
% KWW_Time = KWW_Time_Function(t);
% KWW_Time = gpuArray(KWW_Time);


% s = gpuArray(s);
% 
% LTime = zeros(length(s));
% LKWW = zeros(length(s));
% 
% 
% for i=1:numel(s)
%    LTime(i)=trapz(t,Time.*exp(-s(i)*t)); 
% end
% 
% 
% for i=1:numel(s)
%    LKWW(i)=trapz(t,KWW_Time.*exp(-s(i)*t)); 
% end
% 
% 
% 
% figure()
% hold on
% set(gca, 'XScale', 'log')
% 
% plot(s, LTime);
% plot(s, LKWW);
% plot(s, yTime)
% plot(s, yKWW)
% 
% legend('LTime', 'LKWW', 'yTime', 'yKWW')
% 
% title('Time to Frequency');
% xlabel('log(s)');
% ylabel('|Z|');
