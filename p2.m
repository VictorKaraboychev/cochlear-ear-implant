INPUT_FILE = "input/frequency-sweep.wav";

min_frequency = 100;    % min frequency
max_frequency = 8000;   % max frequency
num_buckets = 4;        % number of buckets

order = 500;            % filter order

error = 0.05;           % frequency error for stopband

bucket_sizes = compute_bucket_sizes(min_frequency, max_frequency, num_buckets);

disp(bucket_sizes);

[audio, sample_rate] = audioread(INPUT_FILE);

% Create a figure for the tile layout
tile_fig = figure;

for i = 1:num_buckets
    f_low = bucket_sizes(i);     
    f_high = bucket_sizes(i + 1);   
    
    filtered = bandpass_filter_fir(audio, f_low, f_high, sample_rate, order);

    % Create subplot in the tile layout for filtered audio (left column)
    subplot(num_buckets, 2, 2*i-1);
    
    % Plot the filtered audio in each bucket
    plot((0:length(filtered)-1)/sample_rate, filtered);
    title(['Filtered ' num2str(i)]);
    xlabel('Time (s)');
    ylabel('Amplitude');

    rectified = abs(filtered);

    amplitude = lowpass_filter_fir(rectified, 400, sample_rate, order);
    
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

function bucket_sizes = compute_bucket_sizes(f_min, f_max, num_buckets)
    linearModel = @(x) f_min + x * (f_max - f_min);
    sqrtModel = @(x, n) f_min + x.^n * (f_max - f_min);
    exponentialModel = @(x) f_min * exp(x * log(f_max/f_min));

    bucket_sizes = exponentialModel(linspace(0, 1, num_buckets + 1));
    disp(bucket_sizes);
end

function filtered_audio = bandpass_filter_irr(audio, f_low, f_high, f_sample, error)
    % f_low: First Passband Frequency
    % f_high: Second Passband Frequency

    % f_sample: Sampling Frequency
    
    tol = 0.0;

    order = 100;                                % Order
    stop_low = f_low * (1 - error + tol);             % First Stopband Frequency
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