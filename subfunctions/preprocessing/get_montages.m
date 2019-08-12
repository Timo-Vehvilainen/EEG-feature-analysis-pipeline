%Timo Vehviläinen August 2019

function [predrug_montages, postdrug_montages, full_montages, obj_array] ...
    = get_montages(data_path, drug_string)
%GET_MONTAGES Imports data from Nicolet .e-files into channel-montaged
%arrays, divided into predrug and postdrug sections
%   PARAMETERS:
%       data_path:
%           a string, representing the full path to the folder containing
%           the .e-files to be imported
%       drug_string:
%           a string representing the name of the drug being analysed. The
%           first appearance of this string is searched for in the
%           annotations of the .e-files. Samples before the annotated moment 
%           are sectioned as "predrug" and samples after it are sectioned
%           as "postdrug".
%
%OUTPUT:
%   All outputs (except sampling_rate) are cell vectors, where each cell
%   element represents the montage data for one subject/patient.
%   All the output montages are 2D-arrays of the form [time_samples * channels],
%   and all of the montages contain the following 8 columns (4 monopolar &
%   4 bipolar channels) in order:
%       [F3, F4, P3, P4, F3-P3 (Left), F4-P4 (Right), F3-F4 (Frontal), P3-P4 (Parietal)]   
%           predrug_montages:
%               contains the channel-montaged samples from before the first
%               appearance of the drug_string - annotation in the .e-files  
%           postdrug_montages:
%               contains the channel-montaged samples from after the first
%               appearance of the drug_string - annotation in the .e-files
%           full_montages:
%               a combination of the predrug_montages and
%               postdrug_montages, containing the full channel-montaged
%               versions of the EEG data in the .e-files.
%           obj_array:
%               a cell vector, where each element is a Nicolet-object
%               struct, which contains all the available data for the given
%               .e-file

data_files = dir(data_path);
all_subs = numel(data_files);

[full_montages, predrug_montages, postdrug_montages, ...
    obj_array] = deal(cell(all_subs, 1));

sampling_rate = 250;

f = waitbar(0, 'Importing EEG data...');
for sub = 1:all_subs
    waitbar(sub/all_subs, f, sprintf('Importing EEG data (subject %d/%d)...', sub, all_subs));
    filename = strcat('Sofia_fentanyl_data/', data_files(sub).name);

    try
        obj = NicoletFile(filename);
    
    %if there were errors while reading the file, skip that patient
    catch e %e is an MException struct
        fprintf(1, 'ERROR! \n');
    end
    
    obj_array{sub} = obj;
    segment_no = numel(obj.segments);
    
    eventMarkers = obj.eventMarkers;
    
    %Assume that the sampling rate is the same for all 4 channels, using F3
    %as the reference.
    annotations = {eventMarkers(:).annotation};
    %replacing empty annotations with 'nothing'
    annotations(find(cellfun(@isempty, annotations))) = {'nothing'};
    
    %Search for the time stamp of the drug administration
    drug_admin = strfind(annotations, drug_string);
    drug_event_index = find(not(cellfun('isempty',drug_admin)));
    if numel(drug_event_index) == 0
        warning('No annotation for "%s" was found in subject no. %d', ...
            drug_string, sub);
        
    end
    drug_event_time = eventMarkers(drug_event_index).dateOLE;
    
    %Flag variable for when we reach the annotation for drug administration
    drug_flag = 0;
    
    full_data = [];
    predrug_data = [];
    postdrug_data = [];
    
    for segment = 1:segment_no

        NrSamples = obj.getNrSamples(segment);
        
        %We are only interested in the F3, F4, P3, P4 channels (the first four)
        e_data = obj.getdata(segment, [1, NrSamples(1)], 1:4);

        %pad any time-gaps from the previous segment with zeros 
        if segment > 1
            time_gap_OLE = obj.segments(segment).dateOLE - ...
                (obj.segments(segment-1).dateOLE + (obj.segments(segment-1).duration / (60*60*24)));
            time_gap_samples = floor(time_gap_OLE * 24 * 60 * 60 * sampling_rate);
        else 
            time_gap_samples = 0;
        end
        
        if time_gap_samples > 0
            full_data =    [full_data; ...
                            zeros(time_gap_samples, 8)];
        end
        
        %Montage data is stored in the form (compired of column vectors):
        % [F3, F4, P3, P4, F3-P3 (Left), F4-P4 (Right), F3-F4 (Frontal), P3-P4 (Posterior)]
        full_data =    [full_data; ...
                        e_data(:, 1:4), ...
                        (e_data(:, 1) - e_data(:, 3)), ...
                        (e_data(:, 2) - e_data(:, 4)), ...
                        (e_data(:, 1) - e_data(:, 2)), ...
                        (e_data(:, 3) - e_data(:, 4))];
        %The variables below hold the same montage data, but separated into
        %predrug and postdrug samples
        
        if drug_flag == 0
            if time_gap_samples > 0
                predrug_data =    [predrug_data; ...
                                    zeros(time_gap_samples, 8)];
            end
            
            %check whether there are still segments to go before drug
            %administration or end of experiment (in case no drug
            %administration annotation was found)
            if (segment < segment_no) && (obj.segments(segment + 1).dateOLE < drug_event_time)
                predrug_data = [predrug_data; ...
                            e_data(:, 1:4), ...
                            (e_data(:, 1) - e_data(:, 3)), ...
                            (e_data(:, 2) - e_data(:, 4)), ...
                            (e_data(:, 1) - e_data(:, 2)), ...
                            (e_data(:, 3) - e_data(:, 4))];
            %else do the last predrug segment, setting the event flag to
            %start filling the postdrug montages
            else
                drug_flag = 1;
                drug_event_sample = round((drug_event_time - obj.segments(segment).dateOLE)*24*60*60*250);
                
                predrug_data = [predrug_data; ...
                            e_data(1:drug_event_sample, 1), e_data(1:drug_event_sample, 2), ...
                            e_data(1:drug_event_sample, 3), e_data(1:drug_event_sample, 4), ...
                            (e_data(1:drug_event_sample, 1) - e_data(1:drug_event_sample, 3)), ...
                            (e_data(1:drug_event_sample, 2) - e_data(1:drug_event_sample, 4)), ...
                            (e_data(1:drug_event_sample, 1) - e_data(1:drug_event_sample, 2)), ...
                            (e_data(1:drug_event_sample, 3) - e_data(1:drug_event_sample, 4))];
                        
                % record any postdrug samples left in the last predrug segment      
                postdrug_data = [postdrug_data;
                        e_data(drug_event_sample+1:end, 1), e_data(drug_event_sample+1:end, 2), ...
                        e_data(drug_event_sample+1:end, 3), e_data(drug_event_sample+1:end, 4), ...
                        (e_data(drug_event_sample+1:end, 1) - e_data(drug_event_sample+1:end, 3)), ...
                        (e_data(drug_event_sample+1:end, 2) - e_data(drug_event_sample+1:end, 4)), ...
                        (e_data(drug_event_sample+1:end, 1) - e_data(drug_event_sample+1:end, 2)), ...
                        (e_data(drug_event_sample+1:end, 3) - e_data(drug_event_sample+1:end, 4))];
            end
        else
            if time_gap_samples > 0
                postdrug_data =    [postdrug_data; ...
                                    zeros(time_gap_samples, 8)];
            end
            
            postdrug_data = [postdrug_data; ...
                            e_data(:, 1:4), ...
                            (e_data(:, 1) - e_data(:, 3)), ...
                            (e_data(:, 2) - e_data(:, 4)), ...
                            (e_data(:, 1) - e_data(:, 2)), ...
                            (e_data(:, 3) - e_data(:, 4))];
        end
    end
    
    full_montages{sub} = full_data;
    predrug_montages{sub} = predrug_data;
    postdrug_montages{sub} = postdrug_data;
end
% fprintf(1, '\nDone.\n');
close(f);
end
