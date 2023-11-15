INPUT_FILE = "input/silent-observers.wav";

min_frequency = 100;    % min frequency
max_frequency = 8000;   % max frequency
num_buckets = 8;        % number of buckets

error = 0.05;           % frequency error for stopband

bucket_sizes = compute_bucket_sizes(min_frequency, max_frequency, num_buckets);

disp(bucket_sizes);

[audio, sample_rate] = audioread(INPUT_FILE);

% Create a figure for the tile layout
tile_fig = figure;

for i = 1:num_buckets
    f_low = bucket_sizes(i);     
    f_high = bucket_sizes(i + 1);   
    
    filtered = bandpass_filter(audio, f_low, f_high, sample_rate, error);

    % Create subplot in the tile layout for filtered audio (left column)
    subplot(num_buckets, 2, 2*i-1);
    
    % Plot the filtered audio in each bucket
    plot((0:length(filtered)-1)/sample_rate, filtered);
    title(['Filtered ' num2str(i)]);
    xlabel('Time (s)');
    ylabel('Amplitude');

    rectified = abs(filtered);

    amplitude = lowpass_filter(rectified, 400, sample_rate, error);
    
    % Create subplot in the tile layout for lowpass audio (right column)
    subplot(num_buckets, 2, 2*i);
    
    % Plot the lowpass audio in each bucket
    plot((0:length(amplitude)-1)/sample_rate, amplitude);
    title(['Low Pass ' num2str(i)]);
    xlabel('Time (s)');
    ylabel('Amplitude');

    % Save filtered audio
    audiowrite("output/bucket_" + i + ".wav", amplitude, sample_rate);
end

% Adjust layout for better appearance
sgtitle('Audio in Each Bucket');

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
