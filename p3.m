INPUT_FILE = "input/at_least_it_was_here.wav";

REFINE_FILE = "refine/at_least_it_was_here.wav";
OUTPUT_FILE = "output/at_least_it_was_here.wav";

[audio, sample_rate] = process_audio(INPUT_FILE, 16000);
duration = size(audio) / sample_rate;

audiowrite(REFINE_FILE, audio, sample_rate);

min_frequency = 100;    % min frequency
max_frequency = 7500;   % max frequency
num_buckets = 128;      % number of buckets

error = 0.05;           % frequency error for stopband

bucket_sizes = compute_bucket_sizes(min_frequency, max_frequency, num_buckets);

output_audio = zeros(size(audio));

for i = 1:num_buckets
    f_low = bucket_sizes(i);     
    f_high = bucket_sizes(i + 1);   
    
    filtered = bandpass_filter(audio, f_low, f_high, sample_rate, error);
    % filtered = new_bandpass_filter(audio, f_low, f_high, sample_rate);

    rectified = abs(filtered);
    amplitude = lowpass_filter(rectified, 400, sample_rate, error);

    audiowrite("output/bucket_" + i + ".wav", filtered, sample_rate);

    frequency = generate_frequency(sample_rate, duration, sqrt(f_low * f_high));

    output_audio = output_audio + amplitude .* frequency;
end

audiowrite(OUTPUT_FILE, output_audio, sample_rate);

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


function bucket_sizes = compute_bucket_sizes(min_freq, max_freq, num_buckets)
    ratio = max_freq / min_freq;
    factor = ratio^(1 / (num_buckets));
    bucket_sizes = min_freq * factor.^(0:(num_buckets));
end

function filtered_audio = bandpass_filter(audio, f_low, f_high, f_sample, error)
    % f_low: First Passband Frequency
    % f_high: Second Passband Frequency

    % f_sample: Sampling Frequency

    stop_low = f_low * (1 - error);             % First Stopband Frequency
    stop_high = f_high * (1 + error);           % Second Stopband Frequency
    attenuation = 150;                          % Stopband Attenuation (dB)
    a_pass  = 0.01;                             % Passband Ripple (dB)
    match  = 'passband';                        % Band to match exactly

    h  = fdesign.bandpass(stop_low, f_low, f_high, stop_high, attenuation, a_pass, attenuation, f_sample);
    Hd = design(h, 'ellip', 'MatchExactly', match);

    filtered_audio = filter(Hd, audio);
end

function filtered_audio = new_bandpass_filter(audio, Fc1, Fc2, Fs)
    N = 500;         % Order
    flag = 'scale';  % Sampling Flag

    % Create the window vector for the design algorithm.
    win = blackmanharris(N+1);

    % Calculate the coefficients using the FIR1 function.
    b  = fir1(N, [Fc1 Fc2]/(Fs/2), 'bandpass', win, flag);
    Hd = dfilt.dffir(b);

    filtered_audio = filter(Hd, audio);
end

function filtered_audio = lowpass_filter(audio, f_high, f_sample, error)
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