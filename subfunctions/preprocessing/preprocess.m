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
%           the lowpass threshold to be applied to the data as a 7th order
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

end
