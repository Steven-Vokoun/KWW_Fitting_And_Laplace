num_sims = 2000;
data_size = 10;
mu = 1;

X = makedist('Normal','mu',mu,'sigma',2);

num_miss_true_mean = 0;
bs_averages = zeros(1, data_size, 'gpuArray'); % Preallocate the array
for sim = 1:num_sims
    Xvals = random(X,1,data_size);
    avg = mean(Xvals);

    % Now do bootstrap sampling
    num_bs_sims = 1000;
    for bs_sim = 1:num_bs_sims
        bs_sample = Xvals(randi(data_size,1,data_size));
        bs_avg = mean(bs_sample);
        bs_averages(bs_sim) = bs_avg;
    end

    percentiles = prctile(gather(bs_averages), [2.5, 97.5]);
    if percentiles(1) > mu || percentiles(2) < mu
        num_miss_true_mean = num_miss_true_mean + 1;
    end
end

fprintf('The probability a bootstrap CI does not contain the true mean is %0.2f\n', num_miss_true_mean/num_sims);



% % Define the parameters
% R1 = 10000;
% R2 = 100000 - R1;
% C1 = 22E-9;
% C2 = C1;
% frequencies = logspace(0, 10, 100); 
% 
% % Calculate the equation
% s = 1i * 2 * pi * frequencies;
% numerator = -R2 ./ (s * C2 * R2 + 1);
% denominator = R1 ./ (s * C1 * R1 + 1);
% TF = numerator ./ denominator;
% 
% % Plot the graph
% loglog(frequencies, abs(TF));
% xlabel('Frequency (Hz)');
% ylabel('Magnitude');
% title('Magnitude Plot of the Transfer Function');
% grid on;