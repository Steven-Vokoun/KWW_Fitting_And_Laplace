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



%% Differention of the data




%KWW
N = .5594;
S = 2044;
P = 313420;
KWW_C = 1.6471e-08;

%Real EIS FIT
C = 1.63E-8;
Rp = 2.3E5;
Rs = 2730;

%variables
t =1E-15:1E-7:1;
f = logspace(1,5,500); %10 to 100k
w = 2*pi*f;
w = gpuArray(w);
t = gpuArray(t);
V = .1;

%EIS
Electrode1 = importfile('0(3).txt');
Z_CPE = sqrt((Electrode1.Z1).^2 + (Electrode1.Z2).^2);


%% KWW
KWW_Time_Function = @(t)(((V./S).*exp(-1.*(t./((S.*P.*KWW_C)./(S+P))).^N)) + ((V./(S+P))*(1-(exp(-(t./((S.*P.*KWW_C)./(S+P))).^N)))));
KWW_Time = KWW_Time_Function(t);

Z_KWW = calculate_impedance(KWW_Time, w, t, V);

%% Randles
Time_Function = @(t)(V/Rs)*exp(-t/((Rs*Rp*C)/(Rs+Rp))) + (V/(Rs+Rp))*(1-exp(-t/((Rs*Rp*C)/(Rs+Rp))));
Randles_Time = Time_Function(t);
 
Z_Randles = calculate_impedance(Randles_Time, w, t, V);

[RAW, TT] = German_Test();
figure()


%% Plot
figure()

% First subplot
subplot(1,2,1)
hold on
set(gca, 'XScale', 'log')

% Use different shades of gray for lines
plot(f, abs(Z_KWW), 'Color', [0.3 0.3 0.3], 'LineStyle', '--', 'DisplayName', 'Laplace of the KWW')
plot(f, abs(Z_Randles), 'Color', [0.5 0.5 0.5], 'LineStyle', ':', 'DisplayName', 'Laplace of a Real Capacitor EIS Fit')
plot((Electrode1.Freq), Z_CPE, 'Color', [0.7 0.7 0.7], 'DisplayName', 'RAW EIS Data')
plot(f, abs(RAW), 'Color', [0.2 0.2 0.2], 'DisplayName', 'Laplace of the KWW')

% Title and labels
title('Time to Frequency Bode Plot')
xlabel('Frequency (Hz)')
ylabel('|Z|')

% Set aspect ratio to square
pbaspect([1 1 1]);

% Second subplot
subplot(1,2,2)
hold on

% Use different shades of gray for lines
plot(real(Z_KWW), -imag(Z_KWW), 'Color', [0.3 0.3 0.3], 'LineStyle', '--', 'DisplayName', 'Laplace of the KWW')
plot(real(Z_Randles), -imag(Z_Randles), 'Color', [0.5 0.5 0.5], 'LineStyle', ':', 'DisplayName', 'Laplace of a Real Capacitor EIS Fit')
plot(Electrode1.Z1, Electrode1.Z2, 'Color', [0.7 0.7 0.7], 'DisplayName', 'RAW EIS Data')
plot(real(RAW), -imag(RAW), 'Color', [0.2 0.2 0.2], 'DisplayName', 'Laplace of the KWW')

% Title and labels
title('Time to Frequency Nyquist Plot')
xlabel('Real Impedance')
ylabel('Imaginary Impedance')

% Set aspect ratio to square
pbaspect([1 1 1]);

% Create a single legend with two rows below the subplots
lgd = legend('Location', 'best');
set(lgd, 'Units', 'normalized', 'Position', [0.3, 0.1, 0.4, 0.05], 'Orientation', 'horizontal', 'NumColumns', 2);




function Z = calculate_impedance(Function, w, t, V)
    % Function: The function for which you want to calculate impedance (KWW_Time)
    % Function_Derivative: The derivative of the function (KWW_Prime)
    % w: Array of angular frequencies
    % t: Time vector
    % V: Some constant value
    
    Function_Derivative = gradient(Function);
    % Calculate A_real
    A_real = zeros(size(w));
    for i = 1:numel(w)
        A_real(i) = trapz(Function_Derivative .* cos(w(i) * t));
    end
    A_real = ((1/V)*(A_real)) + (Function(1)/V);
    
    % Calculate A_imag
    A_imag = zeros(size(w));
    for i = 1:numel(w)
        A_imag(i) = trapz(Function_Derivative .* sin(w(i) * t));
    end
    A_imag = (-1/V) * A_imag;
    
    A = complex(A_real, A_imag);
    
    % Create the complex impedance vector Z
    Z = 1./A;
end



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



function y_dot = robustDiff(y, dt, N)
%ROBUSTDIFF differentiates using smooth noise-robust differentiation formula
%
%   y_dot = robustDiff(y, dt, N)
%
%% Inputs
% y  - signal/vector to differentiate
% dt - time or distance between points in y
% N  - Number of points to use in differentiation.  This value must be
%      positive odd integer greater than or equal 5.
%
%% Outputs
% y_dot - differentiated signal/vector
%
%% Description
% robustDiff differentiates a signal/vector numerically using N
% points.  Both future information and past information are used to
% calculate the derivative.  In signal processing, this is called non-causal.
% The larger the value of N, the more high frequency noise is suppressed
% unlike Savitsky-Golay filters and Central Difference methods (see references).  
% Note that the derivative is not estimated at the edges of y.  This means that
% (N-1)/2 points at the beginning and end of y_dot are NaN.  See the example.
%
%% Example
%   dt = 0.001; % sampling rate of 1000Hz 
%   t = 0:dt:3; % sec
%   noiseFrequency = 400; % Hz
%   noise = 10*rand(size(t)); % Noise is 10% of signal
%   frequency = 1; %Hz
%   y = 100*sin(2*pi*frequency*t) + noise;
%   N = 21; % Number of points to use to estimate derivative
%   y_dot_estimate = robustDiff(y, dt, N);
%   y_dot_actual = 100*2*pi*frequency*cos(2*pi*frequency*t);
%   subplot(211);
%   plot(t, y);
%   title('y vs. t');
%   subplot(212);
%   plot(t, y_dot_actual, 'DisplayName', 'y''_{actual} of sin(t)'); 
%   hold('all')
%   plot(t, y_dot_estimate, 'DisplayName', 'y''_{estimate} of sin(t) + noise');
%   legend('show');
%   hold('off');
%   disp(['Beginning (N-1)/2 points and ending (N-1)/2 points of ' ...
%         'y_dot_estimate are NaN']);
%   y_dot_estimate(1:(N-1)/2)
%   y_dot_estimate(end-(N-1)/2+1:end)
%
%% References
% This function is based on the formulas by Pavel Holoborodko from his
% website: http://www.holoborodko.com/pavel/numerical-methods/numerical-derivative/smooth-low-noise-differentiators/
% A big thanks is due to Pavel Holoborodko for developing these formulas.
%
%   Jason H. Nicholson jashale@yahoo.com
%   $Revision: 1.1 $  $Date: 2014/03/02 16:15:00 $
%% Input Checking
narginchk(3,3);
% N must be odd and greater than 5
if rem(N,2)==1 && N>=5
    % Do nothing
else
    error('N must be postive, odd integer greater than 5.')
end % end if, check N is odd and greater than 5
yLength = length(y);
% y must be have at least N points to estimate the derivative
if yLength >= N
    % Do nothing
else
    error('y must be greater or equal to N.');
end % check yLength is greater or equal N
%% Calculate Coefficients
% Equation for coefficients
% $$\displaystyle {c_k = \frac{1}{2^{2m+1}}\left[{2m\choose m-k+1}-{2m\choose m-k-1}\right]},\quad \displaystyle{m=\frac{N-3}{2}},\quad M=\frac{N-1}{2}$
%
% See reference for more information.
m = (N-3)/2;
M = (N-1)/2;
k = M:-1:1;
% Note that dividing by 2^(2*m+1) should be a bitshift but I don't know how
% to do this in matlab for type double
c_k = (binomialCoefficient(2*m, m-k+1) - binomialCoefficient(2*m, m-k-1))/2^(2*m+1);
%% Calculate coefficients for filter function
b = [c_k 0 -c_k(end:-1:1)];
%% Filter y
% The filter command only computes with past and present information so a
% shift of elements will be required after this step.
y_dot_intermediate = filter(b,1,y(:));
% divide by dt.  Discard first N-1 elements of y_dot_intermediate
y_dot_intermediate = y_dot_intermediate(N:end)/dt;
% shift elements, replace beginning and ending elements, and reshape to the
% same size as y
y_dot = reshape([NaN(M,1); 
                 y_dot_intermediate; 
                 NaN(M,1)],size(y));
end %end function, robustDiff
%% binomialCoefficient
function coefficients = binomialCoefficient(n, k)
% calculates the binomial coeffiecents given k which is a vector and n
% which is a scalar
% preallocated coefficents
coefficients = zeros(size(k));
% find values of k greater than or equal 0 and less than or equal n
index = k >= 0 & k <= n;
% Only calculate coefficients for logical true elements in index
coefficients(index) = arrayfun(@(x) nchoosek(n, x), k(index));
end




