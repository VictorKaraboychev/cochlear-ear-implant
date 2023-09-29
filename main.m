% which -all resample
clc;

INPUT_FILE = "InputSoundtoCochlearImplantProcessor.wav";
OUTPUT_FILE = "OutputSoundtoCochlearImplantProcessor.wav";

[audio, sample_rate] = process_audio(INPUT_FILE, OUTPUT_FILE, 16000);
duration = size(audio) / sample_rate;

% Play the sound file
sound(sound, sample_rate);

% Generate and play 1kHz cos wave as audio
cos_audio = generate_frequency(sample_rate, duration, 1000);
sound(cos_audio, sample_rate);

function cosine_signal = generate_frequency(sample_rate, duration, frequency)
    t = 0:1/sample_rate:duration;
    cosine_signal = cos(2 * pi * frequency * t);

    figure;
    plot(t, cosine_signal);
    title('Cosine Signal Waveform');
    xlabel('Time (seconds)');
    ylabel('Amplitude');
    xlim([0, 2 / frequency]);  % Show two cycles
    grid on;
end

function [audio, sample_rate] = process_audio(input_file, output_file, target_sample_rate)
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
    
    audiowrite(output_file, audio, sample_rate);
end