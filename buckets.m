% Example usage with the Greenwood function
f_min = 100; % Minimum frequency in Hz
f_max = 12800; % Maximum frequency in Hz
N = 7; % Number of bins

% Specify the Greenwood function as the tonotopic model
% greenwoodModel = @(x) (165.4 * (10.^(2.1*x - 0.06)) - 0.11);
linearModel = @(x) f_min + x * (f_max - f_min);
sqrtModel = @(x) f_min + x.^2 * (f_max - f_min);
exponentialModel = @(x) f_min * exp(x * log(f_max/f_min));
% tanhModel = @(x) f_min + 0.5 * (tanh(2*x - 1) + 1) * (f_max - f_min);

% logarithmicModel, exponentialModel, and doublingBinSizeModel are identical
% tanh and greenwood dont seem right...

% figure;
% bins = calc_bins(f_min, f_max, N, greenwoodModel);

bins = calc_bins(f_min, f_max, N, exponentialModel);
bins = calc_bins(f_min, f_max, N, sqrtModel);
% bins = calc_bins(f_min, f_max, N, tanhModel);
bins = calc_bins(f_min, f_max, N, linearModel);

function bins = calc_bins(min_frequency, max_frequency, numBins, tonotopicModel)
    % Calculate the normalized distance along the cochlea for each bin
    x = linspace(0, 1, numBins+1);

    % Calculate the characteristic frequencies using the specified tonotopic model
    bins = tonotopicModel(x);

    % Calculate the bin sizes as the absolute difference between consecutive characteristic frequencies
    sizes = abs(diff([min_frequency, bins, max_frequency]));

    % Display the characteristic frequencies and bin sizes
    format shortG
    disp('Bins (Hz):');
    disp(bins);
    % disp('Bin Sizes (Hz):');
    % disp(binSizes);

    % Plot the characteristic frequencies
    % nexttile;
    % plot(x, bins, 'o-', 'LineWidth', 2);
    % xlabel('Normalized Distance Along Cochlea');
    % ylabel('Characteristic Frequency (Hz)');
    % title('Tonotopic Model: Characteristic Frequencies');
    % grid on;
end