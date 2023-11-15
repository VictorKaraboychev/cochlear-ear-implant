INPUT_FILE = "input/frequency-sweep.wav";

% Example usage:
min_frequency = 100; % replace with your min frequency
max_frequency = 8000; % replace with your max frequency
num_buckets = 8; % replace with your desired number of buckets

error = 0.05;     % replace with your desired frequency error for stopband

bucket_sizes = compute_bucket_sizes(min_frequency, max_frequency, num_buckets);

disp(bucket_sizes);

[audio, sample_rate] = audioread(INPUT_FILE);

% Create a figure for the tile layout
tile_fig = figure;

for i = 1:num_buckets
    f_low = bucket_sizes(i);     
    f_high = bucket_sizes(i + 1);   
    
    filtered_audio = bandpass_filter(audio, f_low, f_high, sample_rate, error);

    filtered_audio = abs(filtered_audio);

    filtered_audio = bandpass_filter(filtered_audio, 1, 400, sample_rate, error);
    
    % Save filtered audio
    audiowrite("output/bucket_" + i + ".wav", filtered_audio, sample_rate);
    
    % Create subplot in the tile layout
    subplot(num_buckets, 1, i);
    
    % Plot the audio in each bucket
    plot((0:length(filtered_audio)-1)/sample_rate, filtered_audio);
    title(['Bucket ' num2str(i)]);
    xlabel('Time (s)');
    ylabel('Amplitude');
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
    attenuation = 50;                           % Stopband Attenuation (dB)
    a_pass  = 0.01;                              % Passband Ripple (dB)
    match  = 'passband';                        % Band to match exactly

    h  = fdesign.bandpass(stop_low, f_low, f_high, stop_high, attenuation, a_pass, attenuation, f_sample);
    Hd = design(h, 'cheby1', 'MatchExactly', match);

    filtered_audio = filter(Hd, audio);
end