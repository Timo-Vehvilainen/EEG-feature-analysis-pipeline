%Timo Vehviläinen August 2019

function [predrug_artifact_masks, postdrug_artifact_masks, full_artifact_masks] ...
    = get_artifact_mask(full_montages, annotation_path, ...
    predrug_durations, postdrug_durations, sampling_rate, global_gate)
%GET_ARTIFACT_MASK Summary of this function goes here
%   This function takes in a montage set of monopolar & bipolar EEG recordings
%   (as given by the function get_montages() ), and produces a mask highlighting all the
%   missing and malformed data, using annotations stored in cell arrays in
%   the folder specified by annotation_folder_path
%   PARAMETERS:
%       full_montages:
%           a cell vector of channel-montaged data, such as the one
%           provided by the get_montages()-function
%       annotation_path:
%           the full file path to the folder containing .mat-formatted
%           annotations. Each .mat-file needs to contain a cell array,
%           where:
%               -the first column contains the text annotations for the
%                   artifacts/periods of missing data, 
%               - the second column specifies
%                   the starting moment of the annotated event (in format
%               hh:mm:ss), and 
%               - the third column specifies the duration of the
%                   annotated event (in format mm:ss). 
%           If no monopolar channel is
%           specified in the text annotation, it is assumed to apply to all
%           channels. However, the channels F3, F4, P3 and P4 can be
%           specified separately. Artifact annotations are detected by the
%           english word "Artifact" or swedish word "Artefakt".
%       predrug_durations:
%           the time-sample durations of the pre-drug sections for each
%           subject, in a cell vector
%       postdrug_durations:
%           the time-sample durations of the postdrug sections for each
%           subject, in a cell vector
%       sampling_rate:
%           the sampling rate of the eeg-recordings
%       global_gate:
%           a value (in microvolts) for automatically magnitude detecting
%           artifacts. Any sample that exceeds this absolute value is
%           marked as an artifact, even if it isn't hand-annotated.
%OUTPUT
%       predrug_artifact_masks:
%           a cell vector, where each element is a 2D binary array representing
%           the predrug artifact data for one subject. They are of the
%           form [time_samples * channels], where 1s represent samples
%           marked as artifacts or missing data, and 0s representing good
%           data.
%       postdrug_artifact_masks:
%           a cell vector, where each element is a 2D binary array representing
%           the postdrug artifact data for one subject. They are of the
%           form [time_samples * channels], where 1s represent samples
%           marked as artifacts or missing data, and 0s representing good
%           data.
%       full_artifact_masks:
%           a concatenation of predrug_artifact_masks and
%           postdrug_artifact_masks


all_subs = numel(full_montages);
[full_artifact_masks, predrug_artifact_masks, postdrug_artifact_masks] ...
    = deal(cell(all_subs, 1));

annotation_files = dir(annotation_path);

% fprintf(1, 'Subject No: ');
f = waitbar(0, 'Constructing artifact mask...');
for sub = 1:all_subs
    
    waitbar(sub/all_subs, f, sprintf('Constructing artifact mask (subject %d/%d)...', sub, all_subs));
%     fprintf(1, '%d... ', sub);

    data = full_montages{sub};
    annotation_filename = fullfile(annotation_files(sub).folder, ...
        annotation_files(sub).name);

    annotation_data = importdata(annotation_filename);
    artifact_mask = logical(zeros(size(data)));
    
    %in get_montages(), time gaps are padded with zeros. These are also to
    %be treated as missing data.
    artifact_mask(data == 0) = 1;
    
    %samples that reach a certain predefined value are also automatically
    %cast as artifacts, even if they are not hand-annotated
    artifact_mask(abs(data) > global_gate) = 1;

    rows = size(annotation_data, 1);
    for row = 1:rows
        
        %detect either English or Swedish artifact annotations
        if ~isempty(strfind(annotation_data{row, 1}, 'Artifact')) || ...
                ~isempty(strfind(annotation_data{row, 1}, 'Artefakt'))
            
            artifact_start = regexp(annotation_data{row, 2}, ':', 'split');
            artifact_duration = regexp(annotation_data{row, 3}, ':', 'split');

            %Convert the artifact timestamps into number of samples from the
            %beginning of the experiment
            artifact_start_samples = max(1, ...
            ((str2num(artifact_start{end-2}) * 60 * 60 + ...        %hours
                str2num(artifact_start{end-1}) * 60 + ...           %minutes
                str2num(artifact_start{end})) * sampling_rate));    %seconds
            
            artifact_duration_samples = ...
                (str2num(artifact_duration{end-1}) * 60 + ...       %minutes
                str2num(artifact_duration{end})) * sampling_rate;   %seconds

            %If there are no monopolar channels specified, apply artifact mask
            %to all channels. Otherwise go through just the flagged
            %channels and mask them (and the bipolar channels they are a part of)
            if isempty(strfind(annotation_data{row, 1}, 'F3')) && ...
                    isempty(strfind(annotation_data{row, 1}, 'F4')) && ...
                    isempty(strfind(annotation_data{row, 1}, 'P3')) && ...
                    isempty(strfind(annotation_data{row, 1}, 'P4'))
                artifact_mask(artifact_start_samples:(artifact_start_samples ...
                    + artifact_duration_samples), :) = 1;
            else
                %F3 artifacts effect F3, Left & Frontal
                if ~isempty(strfind(annotation_data{row, 1}, 'F3'))
                artifact_mask(artifact_start_samples:(artifact_start_samples ...
                    + artifact_duration_samples), [1, 5, 7]) = 1;
                end
                
                %F4 artifacts effect F4, Right & Frontal
                if ~isempty(strfind(annotation_data{row, 1}, 'F4'))
                artifact_mask(artifact_start_samples:(artifact_start_samples ...
                    + artifact_duration_samples), [2, 6, 7]) = 1;
                end
                
                %P3 artifacts effect P3, Left & Parietal
                if ~isempty(strfind(annotation_data{row, 1}, 'P3'))
                artifact_mask(artifact_start_samples:(artifact_start_samples ...
                    + artifact_duration_samples), [3, 5, 8]) = 1;
                end
                
                %P4 artifacts effect P4, Right & Parietal
                if ~isempty(strfind(annotation_data{row, 1}, 'P4'))
                artifact_mask(artifact_start_samples:(artifact_start_samples ...
                    + artifact_duration_samples), [4, 6, 8]) = 1;
                end
            end
        end
    end
    
    full_artifact_masks{sub} = artifact_mask;
    predrug_artifact_masks{sub} = artifact_mask(1:predrug_durations(sub), :);
    postdrug_artifact_masks{sub} = artifact_mask((end-postdrug_durations(sub)+1):end, :);
end
close(f);
% fprintf(1, '\nDone.\n');
end

