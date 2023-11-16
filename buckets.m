% Example usage with the Greenwood function
min_frequency = 100; % Minimum frequency in Hz
max_frequency = 12800; % Maximum frequency in Hz
numBins = 7; % Number of bins

% Specify the Greenwood function as the tonotopic model
greenwoodModel = @(x) 165.4 * (10.^(2.1*x - 0.06)) - 0.11;
linearModel = @(x) min_frequency + x * (max_frequency - min_frequency);
logarithmicModel = @(x) min_frequency * 10.^(x * log10(max_frequency/min_frequency));
sqrtModel = @(x) min_frequency + x.^2 * (max_frequency - min_frequency);
exponentialModel = @(x) min_frequency * exp(x * log(max_frequency/min_frequency));
tanhModel = @(x) min_frequency + 0.5 * (tanh(2*x - 1) + 1) * (max_frequency - min_frequency);
doublingBinSizeModel = @(x) min_frequency * 2.^(x * log2(max_frequency/min_frequency));




% Call the cochlearImplantFilterBank function with the specified tonotopic model
binSizes = cochlearImplantFilterBank(min_frequency, max_frequency, numBins, doublingBinSizeModel);

function binSizes = cochlearImplantFilterBank(min_frequency, max_frequency, numBins, tonotopicModel)
    % Calculate the normalized distance along the cochlea for each bin
    x = linspace(0, 1, numBins+1);

    % Calculate the characteristic frequencies using the specified tonotopic model
    characteristicFrequencies = tonotopicModel(x);

    % Calculate the bin sizes as the absolute difference between consecutive characteristic frequencies
    binSizes = abs(diff([min_frequency, characteristicFrequencies, max_frequency]));

    % Display the characteristic frequencies and bin sizes
    format shortG
    disp('Bins (Hz):');
    disp(characteristicFrequencies);
    % disp('Bin Sizes (Hz):');
    % disp(binSizes);
end