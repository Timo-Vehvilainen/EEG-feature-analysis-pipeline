%Timo Vehviläinen August 2019

function [epoched_data] = epoch_data(data, epoch_scheme, varargin)
%EPOCH_DATA divides 2d data (in the form of [time_samples * channels]  -- that
% is divided into different subject in a cell array -- 
% into epochs of specified duration.
% PARAMETERS:
%           data: A cell vector, where each element represents the data for
%           one EEG subject and is a 2d-array of the form:
%            [time_samples * channels].

%           epoch_scheme: A string that is supposed to have one of the
%               following two values:
%                   'predrug'       OR          'postdrug'
%               The value of this variable determines whether the epoching
%               procedure is initiated from the end of the end of the
%               provided data or the beginning, which in turn affects the
%               position of any remainder samples after epoching (which are
%               discarded). If this value is 'predrug', epoching is started
%               from the end and remainder samples are discarded from the
%               beginning. If the value is 'postdrug', epoching is started
%               from the beginnined and remainder is discarded from the
%               end.
%
%OPTIONAL PARAMETERS:
%           epoch_length_seconds: The length of the epochs in seconds.
%           Defaults to 120 seconds.
%           
%           epoch_overlap_seconds: The amount that consecutive epochs are desired
%           to overlap. Defaults to 60 seconds.
%           
%           max_minutes: The lenght that the data is trimmed down to before
%           epoching. Defaults to 10 minutes.
%
%           sampling_rate: The sampling rate frequency of the data
%           provided. Defaults to 250 Hz.
%
%Output is cell vector of 3d-array-data in the form where each element
%represents the data for a single subject and is a 3D-array of the form:
%           [time_samples, channels, epochs]

%parse input
if ~any(strcmp({'postdrug', 'predrug'}, epoch_scheme))
    error("Incorrect epoching scheme. Please specify 'predrug' or 'postdrug'.\n");
end

p = inputParser;
%optional parameters and their default values
addParameter(p, 'epoch_length_seconds', 120, @isnumeric);
addParameter(p, 'epoch_overlap_seconds', 60, @isnumeric);
addParameter(p, 'max_minutes', 10, @isnumeric);
addParameter(p, 'sampling_rate', 250, @isnumeric);
parse(p, varargin{:});
S = struct2cell(p.Results);
[epoch_length_seconds, epoch_overlap_seconds, max_minutes, sampling_rate] = deal(S{:});

%initialize variables
all_subs = numel(data);
epoched_data = cell(all_subs, 1);
epoch_length_samples = epoch_length_seconds * sampling_rate;
epoch_overlap_samples = epoch_overlap_seconds * sampling_rate;

%Epochs are stored as cell arrays, where each cell is a different subject
%with their own 3d-array of data. The first dimension of the data is time, 
%the 2nd dimension are the channels, and the 3rd
%dimension are the different epochs. 

% fprintf(1, 'Subject No: ');
f = waitbar(0, 'Epoching data...');
for sub = 1:all_subs
%     fprintf(1, '%d... ', sub);
    waitbar(sub/all_subs, f, sprintf('Epoching data (subject %d/%d)...', sub, all_subs));
    %skip "failed to read" - patients
    if isnan(data{sub})
        epoched_data{sub} = NaN;
        continue;
    end
    
    sample_no = size(data{sub}, 1);
    ch_no = size(data{sub}, 2);
    max_samples = min(max_minutes * 60 * sampling_rate, sample_no);
    max_epochs = floor((max_samples - epoch_overlap_samples) / (epoch_length_samples - epoch_overlap_samples));

%     if rem(max_samples - epoch_overlap_samples, epoch_length_samples - epoch_overlap_samples) == 0
%         max_epochs = max_epochs - 1;
%     end
    epoched_data{sub} = zeros(epoch_length_samples, ch_no, max_epochs);
    
    %handle loops separately for each epoch scheme. their difference is in
    %handling remainder samples after dividing into epochs. In predrug, the
    %epochs are read from end towards the beginning, and the 
    %remainder is discarded from the beginning of the run. In postdrug, the
    %data is read from beginning towards the end, and remaining samples
    %are discarded from the end.
    if epoch_scheme == "predrug"
        epoch = 1;
        epoch_end = sample_no;
        epoch_start = epoch_end - epoch_length_samples; 
        while epoch < max_epochs && epoch_start >= 0
            epoched_data{sub}(:, :, epoch) = data{sub}(epoch_start+1:epoch_end, :);
            epoch = epoch + 1;
            epoch_end = epoch_end - (epoch_length_samples - epoch_overlap_samples);
            epoch_start = epoch_end - epoch_length_samples;
        end
        %put the final epochs in correct chronological order
        epoched_data{sub} = flip(epoched_data{sub}, 3);
    
    elseif epoch_scheme == "postdrug"
        epoch = 1;
        epoch_start = 0; 
        epoch_end = epoch_start + epoch_length_samples;
        while epoch < max_epochs && epoch_end <= sample_no
            epoched_data{sub}(:, :, epoch) = data{sub}(epoch_start+1:epoch_end, :);
            epoch = epoch + 1;
            epoch_start = epoch_start + (epoch_length_samples - epoch_overlap_samples);
            epoch_end = epoch_start + epoch_length_samples;
        end   
    end
end
% fprintf(1, '\nDone.\n');
close(f);
end