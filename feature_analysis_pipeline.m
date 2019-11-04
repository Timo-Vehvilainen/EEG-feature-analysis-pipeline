%Timo Vehvilainen November 2019

clc;
clear;
%% FOLDER & FILE PATHS, OTHER VARIABLES (SET THESE ACCORDINGLY):

%Full path to all of the eeg-files (ending in '*.e')
data_path = fullfile(pwd, 'data', '*.e');

%Full path to the artifact annotation files (ending in '*.mat'),
%See documentation for get_artifact_mask() on further info on
%the format of these files
annotation_path = fullfile(pwd, 'artifact_annotations', '*.mat');

%Full path to the folder where you wish to save the resulting plots
results_path = fullfile(pwd, 'results');

%Create the above-defined result folder if it doesn't exist yet
if ~exist(results_path, 'dir')
    mkdir(results_path)
end

%Set variables for pipeline:

%The string to be searched for the in the annotations of the .e-files,
%corresponding to the drug andministration events. In case there are typos
%in the annotations, it might be a good idea to look into the annotations
%inside obj_array.eventMarkers.annotations, and play with this variable
drug_string = 'drug';

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

%(It is normal to see a warning about multiple TSinfo packets when running this)
fprintf(1, 'Reading data from Nicolet .e-files...\n');
[predrug_montages, postdrug_montages, full_montages, obj_array] ...
    = get_montages(data_path, drug_string);
fprintf(1, 'Done importing data.\n')

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
fprintf(1, 'Done.\n');

fprintf(1, 'Applying preprocessing filters to postdrug data...\n');
processed_postdrug_montages = preprocess(postdrug_montages, ...
    'global_highpass', global_highpass, 'global_lowpass', global_lowpass, ...
    'sampling_rate', sampling_rate, 'global_gate', global_gate);
fprintf(1, 'Done.\n');

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
fprintf(1, 'Done.\n');

fprintf(1, 'Epoching predrug artifact_mask...\n');
predrug_artifact_epochs = epoch_data(predrug_artifact_masks, 'predrug', ...
    'epoch_length_seconds', epoch_length_seconds, ...
    'epoch_overlap_seconds', epoch_overlap_seconds, ...
    'sampling_rate', sampling_rate, 'max_minutes', max_predrug_minutes);
fprintf(1, 'Done.\n');

fprintf(1, 'Epoching postdrug montages...\n');
processed_postdrug_epochs = epoch_data(processed_postdrug_montages, 'postdrug', ...
    'epoch_length_seconds', epoch_length_seconds, ...
    'epoch_overlap_seconds', epoch_overlap_seconds, ...
    'sampling_rate', sampling_rate, 'max_minutes', max_postdrug_minutes);
fprintf(1, 'Done.\n');

fprintf(1, 'Epoching postdrug artifact_mask...\n');
postdrug_artifact_epochs = epoch_data(postdrug_artifact_masks, 'postdrug', ...
    'epoch_length_seconds', epoch_length_seconds, ...
    'epoch_overlap_seconds', epoch_overlap_seconds, ...
    'sampling_rate', sampling_rate, 'max_minutes', max_postdrug_minutes);
fprintf(1, 'Done.\n');
                    
%% Construct feature vectors (OPTION TO SET ARTIFACT PERCENTAGE THRESHOLDS)

%The feature structs are saved into .mat-files. This way after being
%constructed once, they can simply be loaded for time-efficiency. 
%You can switch between running the calculation and loading
%already calculated features by commenting (Ctrl+R) and uncommenting (Ctrl+T) 
%parts of this section.

%The artifact percentage thresholds used to reject ecpohs that have an
%excessive amount of unusable sample data. These are set individually for
%each of the features, using a number in the range [0-1]. 
%   -A 0 means that the feature has zero tolerance for artifacts, 
%       and any epoch containing even a single sample that has been 
%       flagged as an artifact will be omitted.
%   - A value of 1 means that the epoch is omitted only if 100% of its
%       samples are marked as artifacts
%For cross-channel features like ASI and cPSD, the percentage is calculated
%as a temporal union of both channels
artifact_percentage_thresholds = ...
    struct( 'ASI',      0.2, ...
            'wPLI',     0.1, ...
            'NC',       0.1, ...
            'aEEG',     0.5, ...
            'rEEG',     0.5, ...
            'PSD',      0.5, ...
            'cPSD',     0.5, ...
            'MFDFA',    0.01, ...
            'SC',       0.5);

warning on verbose
warning on backtrace
fprintf(1, 'Calculating features for predrug epochs...\n');
predrug_features = get_features(...
    processed_predrug_epochs, ...
    predrug_artifact_epochs, ...
    artifact_percentage_thresholds, ...
    'epoch_length_seconds', epoch_length_seconds, ...
    'sampling_rate', sampling_rate);
 save("predrug_features.mat", 'predrug_features')
 fprintf(1, '\nDone.\n');
 
fprintf(1, 'Calculating features for postdrug epochs...\n');
postdrug_features = get_features(...
    processed_postdrug_epochs, ...
    postdrug_artifact_epochs, ...
    artifact_percentage_thresholds, ...
    'epoch_length_seconds', epoch_length_seconds, ...
    'sampling_rate', sampling_rate);
 save("postdrug_features.mat", 'postdrug_features')
 fprintf(1, '\nDone.\n');

% load('predrug_features.mat');
% load('postdrug_features.mat');

%% Color-coding the percentage of artifacts in epochs per subject (optional)

%Additional visualization of how the artifacts effect feature calculation
plot_artifact_percentages(predrug_artifact_epochs,...
    postdrug_artifact_epochs,...
    artifact_percentage_thresholds,...
    results_path)

%% Plot features (SET THESE VARIABLES ACCORDINGLY)

%Different method options for calculating the baseline in the 
%correlation-delta and baseline timeline plots.
%Switch out for the appropriate one, or construct a new one 
%(needs to take a single numeric vector as an argument):
baseline_fcn = @baseline_median5;   %Median of the last 5 predrug epochs.
% baseline_fcn = @baseline_median3;     %Median of the last 3 predrug_epochs.
% baseline_fcn = @baseline_last_epoch;      %The last predrug epoch as baseline.

%Choose whether the baseline-plots are shown with absolute feature values
%or percentile values relative to the baseline:
baseline_mode = 'relative';
% baseline_mode = 'absolute';

%Maximum amounts of pre- & postdrug epochs that are shown on the plots.
%A vector with two elements. The first element determines the amount of
%predrug epochs plotted, and the second determines the number of postdrug
%epochs.
show_epochs = [5, 10];

%To plot just some of the features, comment others out from the following
%command
features_to_plot = {
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

%Drug abbreviation written in labels & titles in the plots:
drug_abbr = 'FE';

%p-value significance level alpha (correlations with p-value lower than this will
%be highlighted in the plots)
p_threshold = 0.05;

%The type of correlation coefficient to be calculated. Can be one of:
%'Spearman', 'Pearson' or 'Kendall'
corr_type = 'Spearman';

%if false, signranks and correlations will be plotted in a two-sided graph,
%with the p-values showing on the right side axis. If true, only the
%signrank medians and correlations coefficients will show, but those with a
%p-value below the significance level will still be highlighted.
hide_p_values = true;

if ~exist(results_path, 'dir')
    mkdir(results_path)
end

fprintf(1, 'Plotting features...\n');
f = waitbar(0, 'Plotting features...',  'Name','Feature plot progress');

%Loop trough chosen features
for feature_idx = 1:numel(features_to_plot)
    feature = features_to_plot{feature_idx};
    waitbar(feature_idx/numel(features_to_plot), f, strcat('Plotting features... ',feature), ...
        'Name','Feature calculation progress');
    set(f, 'HandleVisibility', 'off');

    plotting_toolbox( ...
        predrug_features, postdrug_features, ...% The feature matrices
        feature, ...                            % Chosen feature to be plotted
        results_path, ...                       % Where to save the plots
        ...                                 % OPTIONAL VARIABLES:
        'baseline_fcn', baseline_fcn, ...   % Function to calculate baseline (Defaults to @baseline_median5)
        'baseline_mode', baseline_mode, ... % plot absolute or relative baseline timeline (defaults to relative)
        'show_epochs', show_epochs, ...     % How many epochs to show in the plots (defaults to [5 10])
        'drug_abbr', drug_abbr, ...         % drug abbreviation used in labels & titles (defaults to 'DRUG')
        'p_threshold', p_threshold, ...     % p-value significance threshold (defaults to 0.05)
        'hide_p_values', hide_p_values, ... % whether to only highlight significant p-values (defaults to false)
        'corr_type', corr_type);            % Type of correlation to calculate (defaults to 'Spearman')
    
    set(f, 'HandleVisibility', 'on');
end
fprintf(1, 'Done.\n');
close(f);
close all;
