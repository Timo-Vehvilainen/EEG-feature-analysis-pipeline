%Timo Vehvilainen August 2019

clc;
clear;


%% FOLDER & FILE PATHS, OTHER VARIABLES (SET THESE THESE ACCORDINGLY):

%Full path to all of the eeg-files (ending in '*.e')
data_path = "\*.e";

%Full path to the aftifact annotation files (ending in '*.mat'),
%See documentation for get_artifact_mask() on further info on
%the format of these files
annotation_path = "\*.mat";

%Full path to the folder where you wish to save the resulting plots
results_path = "\results";

%Set variables for pipeline

%Length of a single epoch, in seconds
epoch_length_seconds = 120;

%How many seconds do consecutive epochs overlap with each other
epoch_overlap_seconds = 60;

%How many minutes of data (both predrug & postdrug) to preserve after
%epoching (at maximum) for feature calculation.
max_predrug_minutes = 20;
max_postdrug_minutes = 20;

%Global preprocessing filter settings:
global_highpass = 0.2;  %[Hz]
global_lowpass = 30;    %[Hz]
global_gate = 1000;     %[microVolts]
sampling_rate = 250;    %[Hz]

%% Import data, extract bipolar montages

%The string to be searched for the in the annotations of the .e-files,
%corresponding to the drug andministration events. In case there are typos
%in the annotations, it might be a good idea to look into the annotations
%inside the obj_array, and play around with this variable:
drug_string = '';

fprintf(1, 'Reading data from Nicolet .e-files...\n');
[predrug_montages, postdrug_montages, full_montages, obj_array] ...
    = get_montages(data_path, drug_string);

%% Costruct artifact masks from the annotations

%get sample durations of predrug & postdrug montages
[predrug_durations, ~] = cellfun(@size, predrug_montages);
[postdrug_durations, ~] = cellfun(@size, postdrug_montages);

fprintf(1, 'Constructing artifact masks...\n');
[predrug_artifact_masks, postdrug_artifact_masks, full_artifact_masks] = ...
        get_artifact_mask(full_montages, annotation_path, ...
        predrug_durations, postdrug_durations, sampling_rate, global_gate);

%% Preprocessing - Global Low pass, High pass, Automated artifact flagging

fprintf(1, 'Applying preprocessing filters to predrug data...\n');
processed_predrug_montages = preprocess(predrug_montages, ...
    'global_highpass', global_highpass, 'global_lowpass', global_lowpass, ...
    'sampling_rate', sampling_rate, 'global_gate', global_gate);

fprintf(1, 'Applying preprocessing filters to postdrug data...\n');
processed_postdrug_montages = preprocess(postdrug_montages, ...
    'global_highpass', global_highpass, 'global_lowpass', global_lowpass, ...
    'sampling_rate', sampling_rate, 'global_gate', global_gate);

%% Chop up predrug and postdrug montages to epochs of standard length

max_predrug_epochs = min(floor((max_predrug_minutes*60) / epoch_length_seconds), ...
    (max_predrug_minutes*60 - (epoch_length_seconds - epoch_overlap_seconds)) / epoch_overlap_seconds);

max_postdrug_epochs = min(floor((max_postdrug_minutes*60) / epoch_length_seconds), ...
    (max_postdrug_minutes*60 - (epoch_length_seconds - epoch_overlap_seconds)) / epoch_overlap_seconds);

fprintf(1, 'Epoching predrug montages...\n');
processed_predrug_epochs = epoch_data(processed_predrug_montages, 'predrug', ...
    'epoch_length_seconds', epoch_length_seconds, ...
    'epoch_overlap_seconds', epoch_overlap_seconds, ...
    'sampling_rate', sampling_rate, 'max_minutes', max_predrug_minutes);

fprintf(1, 'Epoching predrug artifact_mask...\n');
predrug_artifact_epochs = epoch_data(predrug_artifact_masks, 'predrug', ...
    'epoch_length_seconds', epoch_length_seconds, ...
    'epoch_overlap_seconds', epoch_overlap_seconds, ...
    'sampling_rate', sampling_rate, 'max_minutes', max_predrug_minutes);

fprintf(1, 'Epoching postdrug montages...\n');
processed_postdrug_epochs = epoch_data(processed_postdrug_montages, 'postdrug', ...
    'epoch_length_seconds', epoch_length_seconds, ...
    'epoch_overlap_seconds', epoch_overlap_seconds, ...
    'sampling_rate', sampling_rate, 'max_minutes', max_postdrug_minutes);

fprintf(1, 'Epoching postdrug artifact_mask...\n');
postdrug_artifact_epochs = epoch_data(postdrug_artifact_masks, 'postdrug', ...
    'epoch_length_seconds', epoch_length_seconds, ...
    'epoch_overlap_seconds', epoch_overlap_seconds, ...
    'sampling_rate', sampling_rate, 'max_minutes', max_postdrug_minutes);

%% Color-coding the percentage of artifacts in epochs per subject

%Additional visualization of how the artifacts effect feature calculation
plot_artifact_percentages(predrug_artifact_epochs,postdrug_artifact_epochs, results_path)
                    
%% Construct feature vectors

%The feature structs are saved into .mat-files, so that after being
%constructed once, they can simply be loaded. 
%You can switch between running the calculation and loading
%already calculated features by commenting (Ctrl+R) and uncommenting (Ctrl+T) 
%parts of this section.

fprintf(1, 'Calculating features for predrug epochs...\n');
predrug_features = get_features(processed_predrug_epochs, predrug_artifact_epochs);
 save("predrug_features.mat", 'predrug_features')
 fprintf(1, '\nDone.\n');
 
fprintf(1, 'Calculating features for postdrug epochs...\n');
postdrug_features = get_features(processed_postdrug_epochs, postdrug_artifact_epochs);
 save("postdrug_features.mat", 'postdrug_features')
 fprintf(1, '\nDone.\n');

% load('predrug_features.mat');
% load('postdrug_features.mat');
%% Plot features

%Different method options for calculating the baseline. Switch out for the
%appropriate one, or construct a new one (needs to take a single numeric
%vector as an argument):

baseline_fcn = @baseline_median5;   %Median of the last 5 predrug epochs.
% baseline_fcn = @baseline_median3;     %Median of the last 3 predrug_epochs.
% baseline_fcn = @baseline_last_epoch;      %The last predrug epoch as baseline.

%Choose whether the baseline-plots are shown with absolute feature values
%or percentile values relative to the baseline:

baseline_mode = 'relative';
% baseline_mode = 'absolute';

%Maximum amounts of pre- & postdrug epochs that are shown on the plots.
show_epochs = 10;

%To plot just some of the features, comment others out from the following
%command
features = {
    'ASI';
    'wPLI';
    'NC';
    'PSD';
    'cPSD';
    'aEEG';
    'rEEG';
    'SC';
    'MFDFA'
};

%Drug abbreviation used in plots:
drug_abbr = 'DRUG';

%p-value significance level (correlations with p-value lower than this will
%be highlighted)
p_threshold = 0.05;

if ~exist(results_path, 'dir')
    mkdir(results_path)
end

fprintf(1, 'Plotting features...\n');
f = waitbar(0, 'Plotting features...',  'Name','Feature plot progress');
for feature_idx = 1:numel(features)
    feature = features{feature_idx};
%     fprintf('%s, ...\n',feature)
    waitbar(feature_idx/numel(features), f, strcat('Plotting features... ',feature), ...
        'Name','Feature calculation progress');
    set(f, 'HandleVisibility', 'off');

    plotting_toolbox(predrug_features, postdrug_features, ...
        feature, results_path, ...
        'baseline_fcn', baseline_fcn, 'baseline_mode', baseline_mode, ...
        'show_epochs', show_epochs, ...
        'drug_abbr', drug_abbr, ...
        'p_threshold', p_threshold);
    
    set(f, 'HandleVisibility', 'on');
end
fprintf(1, 'Done.\n');
close(f);
close all;
