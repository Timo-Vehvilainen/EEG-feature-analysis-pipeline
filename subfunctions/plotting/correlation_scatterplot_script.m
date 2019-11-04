%% Plotting a single correlation scatter plot

%This script is separate from the main pipeline, and
%provides additional functionality for plotting the correlation
%scatter plot of any single postdrug epoch of any feature. The selection
%process is handled with dialog boxes. Before running this, you need to have
%predrug_features and postdrug_features and drug_abbr in your workspace.


baseline_fcn = @baseline_median5;   %Median of the last 5 predrug epochs.
% baseline_fcn = @baseline_median3;     %Median of the last 3 predrug_epochs.
% baseline_fcn = @baseline_last_epoch;      %The last predrug epoch as baseline.

feature_names = fieldnames(predrug_features);
channel_names = {'F3', 'F4', 'P3', 'P4', 'Left', 'Right', 'Frontal', 'Parietal'};
[feature_index, ~] = listdlg('PromptString', 'Select a feature',...
    'ListString',feature_names, 'SelectionMode','single');

feature = feature_names{feature_index};
PSD_band_names = {'\delta-waves (1-3 Hz)', ...
    '\theta-waves (3-8 Hz)', ...
    '\alpha-waves (8-15 Hz)', ...
    '\beta-waves (15-30 Hz)'};
NC_band_names = {'\theta-waves (3-8 Hz)', ...
    '\alpha-waves (8-15 Hz)', ...
    '\beta-waves (15-30 Hz)'};
wPLI_band_names = {'\delta-waves (0.4-3 Hz)',...
    '\theta-waves (3-8 Hz)', ...
    '\alpha-waves (8-15 Hz)', ...
    '\beta-waves (15-30 Hz)'};

[pre, post] = deal(cell(numel(predrug_features), 1));
band_name = '';
ch_name = '';

if strcmp(feature, 'ASI')
    [ch1, ~] = listdlg('Promptstring', 'Select a channel',...
        'ListString',channel_names, 'SelectionMode','single');
    [ch2, ~] = listdlg('Promptstring', 'Select a second channel',...
    'ListString',channel_names, 'SelectionMode','single');
    ch_name = sprintf(', %s vs %s', channel_names{ch1}, channel_names{ch2});
    for sub = 1:numel(predrug_features)
        pre{sub} = predrug_features(sub).(feature)(ch1, ch2, :);
        post{sub} = postdrug_features(sub).(feature)(ch1, ch2, :);
    end
elseif strcmp(feature, 'SC')
    for sub = 1:numel(predrug_features)
        pre{sub} = predrug_features(sub).(feature)';
        post{sub} = postdrug_features(sub).(feature)';
    end
elseif strcmp(feature, 'cPSD')
    [ch1, ~] = listdlg('Promptstring', 'Select a channel',...
        'ListString',channel_names, 'SelectionMode','single');
    [ch2, ~] = listdlg('Promptstring', 'Select a second channel',...
    'ListString',channel_names, 'SelectionMode','single');
    [band, ~] = listdlg('Promptstring', 'Select a frequency band',...
    'ListString',PSD_band_names, 'SelectionMode','single');
    band_name = sprintf(', %s', PSD_band_names{band});
    ch_name = sprintf(', %s vs %s', channel_names{ch1}, channel_names{ch2});
    for sub = 1:numel(predrug_features)
        pre{sub} = predrug_features(sub).(feature)(ch1, ch2, band, :);
        post{sub} = postdrug_features(sub).(feature)(ch1, ch2, band, :);
    end
elseif strcmp(feature, 'wPLI')
    [ch1, ~] = listdlg('Promptstring', 'Select a channel',...
        'ListString',channel_names, 'SelectionMode','single');
    [ch2, ~] = listdlg('Promptstring', 'Select a second channel',...
    'ListString',channel_names, 'SelectionMode','single');
    [band, ~] = listdlg('Promptstring', 'Select a frequency band',...
    'ListString',wPLI_band_names, 'SelectionMode','single');
    band_name = sprintf(', %s', wPLI_band_names{band});
    ch_name = sprintf(', %s vs %s', channel_names{ch1}, channel_names{ch2});
    for sub = 1:numel(predrug_features)
        pre{sub} = predrug_features(sub).(feature)(ch1, ch2, band, :);
        post{sub} = postdrug_features(sub).(feature)(ch1, ch2, band, :);
    end
elseif strcmp(feature, 'NC') 
    [ch1, ~] = listdlg('Promptstring', 'Select a channel',...
        'ListString',channel_names, 'SelectionMode','single');
    [band, ~] = listdlg('Promptstring', 'Select a frequency band',...
    'ListString',NC_band_names, 'SelectionMode','single');
    band_name = sprintf(', %s', NC_band_names{band});
    ch_name = sprintf(', %s', channel_names{ch1});
    for sub = 1:numel(predrug_features)
        pre{sub} = predrug_features(sub).(feature)(ch1, band, :);
        post{sub} = postdrug_features(sub).(feature)(ch1, band, :);
    end
elseif strcmp(feature, 'PSD')
    [ch1, ~] = listdlg('Promptstring', 'Select a channel',...
        'ListString',channel_names, 'SelectionMode','single');
    [band, ~] = listdlg('Promptstring', 'Select a frequency band',...
    'ListString',PSD_band_names, 'SelectionMode','single');
    band_name = sprintf(', %s', PSD_band_names{band});
    ch_name = sprintf(', %s', channel_names{ch1});
    for sub = 1:numel(predrug_features)
        pre{sub} = predrug_features(sub).(feature)(ch1, band, :);
        post{sub} = postdrug_features(sub).(feature)(ch1, band, :);
    end
elseif any(strcmp(feature, {'aEEG', 'rEEG', 'MFDFA'}))
    subfeatures = fieldnames(predrug_features(1).(feature));
    [subfeature_index, ~] = listdlg('Promptstring', 'Select a subfeature',...
        'ListString',subfeatures, 'SelectionMode','single');
    subfeature = subfeatures{subfeature_index};
    [ch1, ~] = listdlg('Promptstring', 'Select a channel',...
        'ListString',channel_names, 'SelectionMode','single');
    band_name = sprintf(', %s', subfeature);
    ch_name = sprintf(', %s', channel_names{ch1});
    for sub = 1:numel(predrug_features)
        pre{sub} = predrug_features(sub).(feature).(subfeature)(ch1, :)';
        post{sub} = postdrug_features(sub).(feature).(subfeature)(ch1, :)';
    end
end

post_epochs = numel(post{1});
postdrug_epoch_cell = inputdlg(...
    sprintf('Give postdrug epoch to calculate correlation (integer in range [1-%d])', post_epochs));
postdrug_epoch = str2double(postdrug_epoch_cell{1});

plot_bl_corr_single_epoch(pre, post, baseline_fcn, postdrug_epoch)

title(sprintf('%s%s%s\nCorrelation (postdrug epoch %d)', feature, band_name, ...
    ch_name, postdrug_epoch))
ylabel(sprintf('\\Delta%s', drug_abbr), 'FontWeight','bold')
xlabel('BL', 'FontWeight','bold')