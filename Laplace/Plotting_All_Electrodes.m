% Define parameters for different electrodes
params = struct('N', [.1, 1, 0.407320281], ...
    'S', [1560, 2120, 2005], ...
    'P', [299000, 412000, 808593], ...
    'KWW_C', [1.70e-08, 6.83e-09, 8.90142e-09]);

eisParams = struct('Rs', [1560, 2120, 2005], ...
    'Rp', [2.3E5, 6.17E+05, 159577.5], ...
    'C', [1.63E-8, 6.64E-09, 9.73565E-09]);

% Define file names for different electrodes
fileNames = {'0(3)A.txt', '0(3)B.txt', '0(3)D.txt'};

t = 1E-15:1E-7:1;
f = logspace(1, 5, 500); % 10 to 100k
w = 2 * pi * f;
V = 0.1;
w = gpuArray(w);
t = gpuArray(t);

figure()
for idx = 1:3
    % Load the corresponding file
    currentFile = importfile(fileNames{idx});

    % Extract parameters for the current electrode
    N = params.N(idx);
    S = params.S(idx);
    P = params.P(idx);
    KWW_C = params.KWW_C(idx);
    Rs = eisParams.Rs(idx);
    Rp = eisParams.Rp(idx);
    C = eisParams.C(idx);

    % EIS
    Z_CPE = sqrt((currentFile.Z1) .^ 2 + (currentFile.Z2) .^ 2);

    % Define KWW function
    KWW_Time_Function = @(t) (((V ./ S) .* exp(-1 .* (t ./ ((S .* P .* KWW_C) / (S + P))) .^ N)) + ((V / (S + P)) * (1 - (exp(-(t ./ ((S .* P .* KWW_C) / (S + P))) .^ N)))));
    KWW_Time = KWW_Time_Function(t);

    % Calculate KWW impedance
    Z_KWW = calculate_impedance(KWW_Time, w, t, V);

    % Define Randles function
    Time_Function = @(t) (V / Rs) * exp(-t / ((Rs * Rp * C) / (Rs + Rp))) + (V / (Rs + Rp)) * (1 - exp(-t / ((Rs * Rp * C) / (Rs + Rp))));
    Randles_Time = Time_Function(t);

    % Calculate Randles impedance
    Z_Randles = calculate_impedance(Randles_Time, w, t, V);

    % Subplot setup
    subplot(2, 3, idx)

    % Use different shades of gray for lines
    hold on
    set(gca, 'XScale', 'log')
    plot(f, abs(Z_KWW), 'Color', [0.3 0.3 0.3], 'LineStyle', '--', 'DisplayName', 'Laplace of the KWW')
    plot(f, abs(Z_Randles), 'Color', [0.5 0.5 0.5], 'LineStyle', ':', 'DisplayName', 'Laplace of a Real Capacitor EIS Fit')
    plot(currentFile.Freq, Z_CPE, 'Color', [0.7 0.7 0.7], 'DisplayName', 'EIS Raw Data')

    % Title and labels
    title(char('A' + idx - 1))
    xlabel('Frequency (Hz)')
    ylabel('|Z|')

    % Set aspect ratio to square
    pbaspect([1 1 1]);

    % Subplot setup
    subplot(2, 3, idx+3)

    % Use different shades of gray for lines
    hold on
    plot(real(Z_KWW), -imag(Z_KWW), 'Color', [0.3 0.3 0.3], 'LineStyle', '--', 'DisplayName', 'Laplace of the KWW')
    plot(real(Z_Randles), -imag(Z_Randles), 'Color', [0.5 0.5 0.5], 'LineStyle', ':', 'DisplayName', 'Laplace of a Real Capacitor EIS Fit')
    plot(currentFile.Z1, currentFile.Z2, 'Color', [0.7 0.7 0.7], 'DisplayName', 'EIS Raw Data')

    % Title and labels
    xlabel('Real Impedance')
    ylabel('Imaginary Impedance')

    % Set aspect ratio to square
    pbaspect([1 1 1]);

end

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