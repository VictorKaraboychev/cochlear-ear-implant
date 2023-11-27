FILE = "frequency-sweep.wav";
% FILE = "tektok.wav";

INPUT_FILE = "input/" + FILE;

REFINE_FILE = "refine/" + FILE;

RECONSTRUCTED_FILE = "output/recon-" + FILE;
OUTPUT_FILE = "output/" + FILE;

[audio, sample_rate] = process_audio(INPUT_FILE, 16000);
duration = size(audio) / sample_rate;

audiowrite(REFINE_FILE, audio, sample_rate);

error = 0.05;           % frequency error for stopband

min_frequency = 100;    % min frequency
max_frequency = 8000 * (1 - error);   % max frequency
num_buckets = 8;        % number of buckets

order = 500;            % filter order

bucket_sizes = compute_bucket_sizes(min_frequency, max_frequency, num_buckets);

reconstructed_audio = zeros(size(audio));
output_audio = zeros(size(audio));

for i = 1:num_buckets
    f_low = bucket_sizes(i);     
    f_high = bucket_sizes(i + 1);   
    
    filtered = bandpass_filter_fir(audio, f_low, f_high, sample_rate, order);

    rectified = abs(filtered);
    amplitude = lowpass_filter_fir(rectified, 400, sample_rate, order);

    audiowrite("output/bucket_" + i + ".wav", filtered, sample_rate);

    frequency = generate_frequency(sample_rate, duration, sqrt(f_low * f_high));

    reconstructed_audio = reconstructed_audio + filtered;
    output_audio = output_audio + frequency .* amplitude;
end

audiowrite(RECONSTRUCTED_FILE, reconstructed_audio, sample_rate);
audiowrite(OUTPUT_FILE, output_audio, sample_rate);

disp("DONE");

function [audio, sample_rate] = process_audio(input_file, target_sample_rate)
    [raw_audio, sample_rate] = audioread(input_file);
    
    shape = size(raw_audio);
    
    samples = shape(1);
    channels = shape(2);
    
    audio = zeros(samples, 1);
    
    % Flattening audio channels to mono
    if channels > 1
        for s = 1:samples
            for c = 1:channels
                audio(s) = audio(s) + raw_audio(s, c);
            end
        end
    end
    
    % Resample the signal to target_sample_rate
    audio = resample(audio, target_sample_rate, sample_rate);
    sample_rate = target_sample_rate;
end

function cosine_signal = generate_frequency(sample_rate, duration, frequency)
    t = (0:1/sample_rate:duration - 1/sample_rate).';
    cosine_signal = cos(2 * pi * frequency * t);
end


function bucket_sizes = compute_bucket_sizes(f_min, f_max, num_buckets)
    linearModel = @(x) f_min + x * (f_max - f_min);
    sqrtModel = @(x, n) f_min + x.^n * (f_max - f_min);
    exponentialModel = @(x) f_min * exp(x * log(f_max/f_min));
    % sinModel = @(x, n) f_min + (f_max - f_min) .* sin((pi/2)*x.^(1/n));
    % sinModel2 = @(x, n) f_min + (f_max - f_min) .* (x - (sin(2*pi.*x + pi) / (2*pi)));
    bucket_sizes = linearModel(linspace(0, 1, num_buckets + 1));
    disp(bucket_sizes);
end

function filtered_audio = bandpass_filter_irr(audio, f_low, f_high, f_sample, error)
    % f_low: First Passband Frequency
    % f_high: Second Passband Frequency

    % f_sample: Sampling Frequency
    
    tol = 0.0;

    stop_low = f_low * (1 - error + tol);       % First Stopband Frequency
    stop_high = f_high * (1 + error);           % Second Stopband Frequency
    attenuation = 50;                           % Stopband Attenuation (dB)
    a_pass  = 0.01;                             % Passband Ripple (dB)
    match  = 'passband';                        % Band to match exactly

    h = fdesign.bandpass(stop_low, f_low * (1 + tol), f_high * (1 - tol), stop_high, attenuation, a_pass, attenuation, f_sample);
    Hd = design(h, 'ellip', 'MatchExactly', match);

    % h = fdesign.bandpass('N,Fp1,Fp2,Ap', order, f_low, f_high, a_pass, f_sample);
    % Hd = design(h, 'cheby1');

    filtered_audio = filter(Hd, audio);
end

function filtered_audio = bandpass_filter_fir(audio, f_low, f_high, f_sample, order)
    flag = 'scale';  % Sampling Flag

    % Create the window vector for the design algorithm.
    win = hamming(order + 1);

    % Calculate the coefficients using the FIR1 function.
    b  = fir1(order, [f_low f_high]/(f_sample/2), 'bandpass', win, flag);
    Hd = dfilt.dffir(b);

    filtered_audio = filter(Hd, audio);
end

function filtered_audio = lowpass_filter_irr(audio, f_high, f_sample, error)
    % f_high: Passband Frequency

    stop = f_high * (1 + error);                % Stopband Frequency
    attenuation = 150;                          % Stopband Attenuation (dB)
    a_pass = 0.01;                              % Passband Ripple (dB)
    match = 'passband';                         % Band to match exactly
    
    % Construct an FDESIGN object and call its CHEBY1 method.
    h  = fdesign.lowpass(f_high, stop, a_pass, attenuation, f_sample);
    Hd = design(h, 'ellip', 'MatchExactly', match);

    filtered_audio = filter(Hd, audio);
end

function filtered_audio = lowpass_filter_fir(audio, f_high, f_sample, order)
    flag = 'scale';                             % Sampling Flag
    
    % Create the window vector for the design algorithm.
    win = hamming(order + 1);
    
    % Calculate the coefficients using the FIR1 function.
    b  = fir1(order, f_high/(f_sample/2), 'low', win, flag);
    Hd = dfilt.dffir(b);

    filtered_audio = filter(Hd, audio);
end