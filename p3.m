FILES = [
    "at_least_it_was_here.wav";
    "frequency-sweep.wav";
    "hotel-california.wav";
    "intrinsically.wav";
    "techies.wav";
    "tektok.wav"
];

BUCKETS = [
    12;
    48;
    128
];

MIN_FREQUENCY = 100;        % min frequency
MAX_FREQUENCY = 8000;       % max frequency

MAX_DB = 0;                 % maximum decibles for normaliztion

PLAY_SOUND = 0;             % play = 1, dont play = 0
for b = 1:length(BUCKETS)
    for f = 1:length(FILES)
        file = FILES(f);
        num_buckets = BUCKETS(b);

        input_file = "input/" + file;
        output_file = "output/" + num_buckets + "buckets_" + file;
        
        % Process the audio
        [audio, sample_rate] = process_audio(input_file, 16000, MIN_FREQUENCY, MAX_FREQUENCY, num_buckets, MAX_DB);
        
        % Play the sound file
        if PLAY_SOUND
            sound(audio, sample_rate);
        end
        
        % Write the audio to a file
        audiowrite(output_file, audio, sample_rate);
    end
end

disp("DONE");

function [audio, sample_rate] = process_audio(input_file, target_sample_rate, min_freq, max_freq, num_buckets, max_db)
    % Read in audio signal
    [raw_audio, sample_rate] = audioread(input_file);
    
    shape = size(raw_audio);
    
    % Get number of samples and channels
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
    duration = size(audio) / sample_rate;
    
    % Generate frequency buckets
    bucket_sizes = compute_bucket_sizes(min_freq, max_freq, num_buckets);
    
    output_audio = zeros(size(audio));
    
    for i = 1:num_buckets
        f_low = bucket_sizes(i);     
        f_high = bucket_sizes(i + 1);   
        
        % Extract frequency band
        filtered = bandpass_filter_rect(audio, f_low, f_high, sample_rate);
    
        % Get amplitude of frequency band
        rectified = abs(filtered);
        amplitude = lowpass_filter_kaiser(rectified, 400, sample_rate);
    
        % Generate frequency and modulate by amplitude
        frequency = generate_frequency(sample_rate, duration, sqrt(f_low * f_high)) .* amplitude;

        % Recombine frequency bands
        output_audio = output_audio + frequency;
    end
    
    % Normalize the output signal
    audio = output_audio * (sqrt(10 .^ (max_db / 10))) / max(abs(output_audio));
    disp("FINISHED PROCESSING FILE: " + input_file);
end

function cosine_signal = generate_frequency(sample_rate, duration, frequency)
    t = (0:1/sample_rate:duration - 1/sample_rate).';
    cosine_signal = cos(2 * pi * frequency * t);
end


function bucket_sizes = compute_bucket_sizes(f_min, f_max, num_buckets)
    linearModel = @(x) f_min + x * (f_max - f_min);
    sqrtModel = @(x, n) f_min + x.^n * (f_max - f_min);
    exponentialModel = @(x) f_min * exp(x * log(f_max/f_min));

    bucket_sizes = exponentialModel(linspace(0, 1, num_buckets + 1));
end

function filtered_audio = bandpass_filter_cheby(audio, f_low, f_high, f_sample)
    % f_low: First Passband Frequency
    % f_high: Second Passband Frequency
    % f_sample: Sampling Frequency
    
    order = 50;
    a_pass  = 0.01;                             % Passband Ripple (dB)

    h = fdesign.bandpass('N,Fp1,Fp2,Ap', order, f_low, f_high, a_pass, f_sample);
    Hd = design(h, 'cheby1');

    filtered_audio = filter(Hd, audio);
end

function filtered_audio = bandpass_filter_kaiser(audio, f_low, f_high, f_sample)
    flag = 'scale';  % Sampling Flag
    order = 500;

    % Create the window vector for the design algorithm.
    win = kaiser(order + 1);

    % Calculate the coefficients using the FIR1 function.
    b  = fir1(order, [f_low f_high]/(f_sample/2), 'bandpass', win, flag);
    Hd = dfilt.dffir(b);

    filtered_audio = filter(Hd, audio);
end

function filtered_audio = bandpass_filter_rect(audio, f_low, f_high, f_sample)
    order = 500;
    flag = 'scale';  % Sampling Flag

    % Create the window vector for the design algorithm.
    win = rectwin(order+1);
    
    % Calculate the coefficients using the FIR1 function.
    b  = fir1(order, [f_low f_high]/(f_sample/2), 'bandpass', win, flag);
    Hd = dfilt.dffir(b);

    filtered_audio = filter(Hd, audio);
end

function filtered_audio = lowpass_filter_cheby(audio, f_high, f_sample)
    % f_high: Passband Frequency
    order = 50;
    a_pass = 0.01;                              % Passband Ripple (dB)
    
    h = fdesign.lowpass('N,Fp,Ap', order, f_high, a_pass, f_sample);
    Hd = design(h, 'cheby1');

    filtered_audio = filter(Hd, audio);
end

function filtered_audio = lowpass_filter_kaiser(audio, f_high, f_sample)
    flag = 'scale';                             % Sampling Flag
    order = 500;
    
    % Create the window vector for the design algorithm.
    win = kaiser(order + 1);
    
    % Calculate the coefficients using the FIR1 function.
    b  = fir1(order, f_high/(f_sample/2), 'low', win, flag);
    Hd = dfilt.dffir(b);

    filtered_audio = filter(Hd, audio);
end

function filtered_audio = lowpass_filter_rect(audio, f_high, f_sample)
    flag = 'scale';  % Sampling Flag
    order = 500;

    % Create the window vector for the design algorithm.
    win = rectwin(order+1);
    
    % Calculate the coefficients using the FIR1 function.
    b  = fir1(order, f_high/(f_sample/2), 'low', win, flag);
    Hd = dfilt.dffir(b);

    filtered_audio = filter(Hd, audio);
end