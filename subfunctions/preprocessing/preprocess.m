function [processed_data_montages] = preprocess(data_montages, varargin)
%PREPROCESS Runs a global band pass filter on the given data montages
%   PARAMETERS:
%       data_montages:
%           a cell vector, where each element represents channel-montaged
%           data for one subject/patient, and is a 2D-array of the form 
%           [time_samples * channels].
%   OPTIONAL PARAMETERS:
%       sampling_rate:
%           the sampling rate frequency of the EEG data. Defaults to 250
%           Hz.
%       global_highpass:
%           the highpass threshold to be applied to the data as a 5th order
%           Butterworth filter (forward-reverse). Defaults to 0.2 Hz.
%       global_lowpass:
%           the highpass threshold to be applied to the data as a 7th order
%           Butterworth filter (forward-reverse). Defaults to 30 Hz.
%       global_gate:
%           a value (in microVolts) for automatically magnitude detecting
%           artifacts. Any sample that exceeds this absolute value is
%           removed and set as 0 before applying the global filters.
%           Defaults to 1000 microVolts

p = inputParser;

%optional variables and their default values
addParameter(p, 'sampling_rate', 250, @isnumeric);
addParameter(p, 'global_highpass', 0.2, @isnumeric);
addParameter(p, 'global_lowpass', 30, @isnumeric);
addParameter(p, 'global_gate', 1000, @isnumeric);

parse(p, varargin{:});
S = struct2cell(p.Results);
[global_gate, global_highpass, global_lowpass, sampling_rate] = deal(S{:});

processed_data_montages = data_montages;
all_subs = numel(data_montages);
ch_no = 8;


f = waitbar(0, 'Preprocessing data...');
for sub = 1:all_subs
    %set all the samples exceeding the global gate to 0. This is to prevent
    %ringing artifact from appearing when filtering
    waitbar(sub/all_subs, f, sprintf('Preprocessing data (subject %d/%d)...', sub, all_subs));
    
    processed_data_montages{sub}(abs(processed_data_montages{sub}) > global_gate) = 0;
   
    for ch = 1:ch_no
        processed_data_montages{sub}(:, ch) = ...
            my_bandpass(processed_data_montages{sub}(:, ch), ...
            [global_highpass, global_lowpass], sampling_rate);
    end
end
close(f);

% all_subs = numel(data);
% preprocessed_data = cell(1, all_subs);
% 
% fprintf(1, 'Subject No: ');
% 
% parfor sub = 1:all_subs
% %     figure(sub);
%     fprintf(1, '%d... ', sub);
%     %skip invalid patients
%     if isnan(data{sub})
%         preprocessed_data{sub} = NaN;
%         continue
%     end
%     
%     active_subject = data{sub};
%     processed_subject = zeros(size(active_subject));
%     
%     channel_names = {'F3', 'F4', 'P3', 'P4', 'Left', 'Right', 'Frontal', 'Parietal'};
%     
% %     ch_no = size(active_subject, 2);
%     ch_no = numel(channel_names);
%     %preprocess each of the montages
%     for ch = 1:ch_no
%         montage_inproc = active_subject(:, ch);
% %         x_axis = 1:numel(montage_inproc);
% %         x_axis = x_axis * (1/sampling_rate);
% %         subplot(3, ch_no, ch)
% %         plot(1:numel(montage_inproc), montage_inproc, 'b');
% %         title(sprintf("Subject %d, channel %s", sub, channel_names{ch}));
%         %montage_inproc = bsxfun(@minus, montage_inproc, mean(montage_inproc)); % remove DC component
%         
%         artifact_mask = artifact_masks{sub}(:, ch);
%         high_amp_artifacts = find(abs(montage_inproc) > global_gate);
%         artifact_mask(high_amp_artifacts) = 1;
%         filtered_parts = montage_inproc;
%         montage_inproc(artifact_mask == 1) = 0;  
%         filtered_parts(artifact_mask == 0) = NaN;
% %         
% %         subplot(3, ch_no, ch_no+ch)
% %         plot(x_axis, montage_inproc, 'b');  
% %         hold on;
% %         plot(x_axis, filtered_parts, 'r');
% %         title("Artifact removal");
% %         xlabel('seconds')
% %         hold off;
% %         
% %         subplot(3, ch_no, 2*ch_no+ch)
% %         plot(x_axis, montage_inproc, 'b');
% %         hold on;
%         
% %         %highpass filter, cutoff is 0.4 Hz
% %         [b,a] = butter(5,highpass/(sampling_rate/2),'high'); 
% %         montage_inproc = filtfilt(b,a,montage_inproc);
% %         plot(x_axis, montage_inproc, 'g');
% %             
% %         %lowpass filter, cutoff is 30 Hz
% %         [b,a] = butter(7,lowpass/(sampling_rate/2),'low'); 
% %         montage_inproc = filtfilt(b,a,montage_inproc);
% %         plot(x_axis, montage_inproc, 'm');
% %         title("Frequency filtering");
% %         legend('Original', 'after high pass', 'after high + low pass')
% %         xlabel('seconds')
% %         hold off;
% 
%         montage_inproc = my_bandpass(montage_inproc, ...
%             [global_highpass, global_lowpass], sampling_rate);
%         
%         processed_subject(:, ch) = montage_inproc;
%     end
%     preprocessed_data{sub} = processed_subject;
%     
% end
% fprintf(1, '\nDone.\n');
end
