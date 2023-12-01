% which -all resample
clc;
tiledlayout(2,1);
INPUT_FILE = "input/audio.wav";
OUTPUT_FILE = "output/audio.wav";

PLAY_SOUND = 1;

[audio, sample_rate] = process_audio(INPUT_FILE, OUTPUT_FILE, 16000);
duration = size(audio) / sample_rate;

% Play the sound file
if PLAY_SOUND
    sound(audio, sample_rate);
end

% plot the sound file
plot_audio(audio, sample_rate);

% Generate and play 1kHz cos wave as audio
cos_audio = generate_frequency(sample_rate, duration, 1000);

if PLAY_SOUND
    sound(cos_audio, sample_rate);
end

function cosine_signal = generate_frequency(sample_rate, duration, frequency)
    t = 0:1/sample_rate:duration;

    cosine_signal = cos(2 * pi * frequency * t);

    nexttile
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


function plot_audio(audio, sample_rate)
    duration = size(audio) / sample_rate;
    
    t = 0:1/sample_rate:duration-(1/sample_rate);

    nexttile
    plot(t, audio);
    title('Audio Waveform');
    xlabel('Time (seconds)');
    ylabel('Amplitude');
    grid on;
end