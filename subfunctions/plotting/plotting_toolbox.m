%Timo Vehviläinen August 2019

function [] = ...
    plotting_toolbox(predrug_features,postdrug_features, feature, results_path, varargin)
%PLOTTING_TOOLBOX This function is used for plotting the data inside a
%single field inside the feature struct produced by the function get_features().
%   Visualizes the feature timelines for all subjects, calculates a
%   baseline with a provided function, visualizes baseline-relative feature
%   values, calculates correlations and saves the results into matlab
%   figures and images.
%   PARAMETERS:
%       predrug_features:
%           a feature struct given be the get_features()-function,
%           representing the feature values for the predrug epochs of all
%           the subjects in the study.
%       postdrug_features:
%           a feature struct given be the get_features()-function,
%           representing the feature values for the postdrug epochs of all
%           the subjects in the study.
%       feature:
%           the feature to be analyzed as a string. Needs to be one of the
%           predefined feature names (specified below in the
%           expected_feature_names variable)
%       results_path:
%           the full file path to the folder where the result figures and
%           images are to be saved in.
%OPTIONAL VARIABLES:
%       epoch_length_seconds: The length of the epochs in seconds.
%           Defaults to 120 seconds.
%       epoch_overlap_seconds: The amount that consecutive epochs are desired
%           to overlap. Defaults to 60 seconds.    
%       show_epochs:
%           a positive number, determining how many pre- and postdrug
%           epochs are visualized in the figures around the point of drug
%           administration. Defaults to 5.
%       sampling_rate: The sampling rate frequency of the data
%           provided. Defaults to 250 Hz.
%       baseline_fcn:
%           a function to be used the calculate the baseline from the
%           predrug epochs. The function needs to take as its parameter a
%           single vector of numeric values. Defaults to @baseline_median5,
%           which calculates the baseline as the median of the last 5
%           predrug epochs.
%       baseline_mode:
%           a string that is either 'relative' or 'absolute'. This
%           parameter determines whether the baseline figures are
%           represented as percentages relative to the baseline (calculated
%           as (postdrug_value / baseline) or as absolute feature value
%           (calculated as postdrug_value - baseline). Defaults to
%           'relative'.
%       drug_abbr:
%           a string representing an abbreviation of the drug being used.
%           It is displayed in the baseline-delta plots. Defaults to
%           "DRUG".
%       p-threshold:
%           The significance level for the correlation calculations. A
%           number in the interval [0, 1]. Defaults to 0.05.

%parse input
expected_channels = {'F3', 'F4', 'P3', 'P4', 'Left', 'Right', 'Frontal', 'Parietal'};
expected_feature_names = {'ASI', 'wPLI', 'NC', 'PSD', 'cPSD', 'aEEG', 'rEEG', 'MFDFA', 'SC'};
if ~any(strcmp(expected_feature_names, feature))
    error(sprintf("Unexpected feature: '%s'. Expected one of:\n", feature) + ...
        sprintf("%s ", expected_feature_names{:}));
end

p = inputParser;
addParameter(p, 'epoch_length_seconds', 120, @isnumeric);
addParameter(p, 'epoch_overlap_seconds', 60, @isnumeric);
addParameter(p, 'show_epochs', 5, @isnumeric);
addParameter(p, 'drug_abbr', 'DRUG', @ischar);
addParameter(p, 'p_threshold', 0.05, @isnumeric);

fcn_validation = @(x)  isa(x, 'function_handle');
addParameter(p, 'baseline_fcn', @baseline_median5, fcn_validation);

expected_baseline_modes = {'relative', 'absolute'};
addParameter(p, 'baseline_mode', 'relative', ...
    @(x) any(validatestring(x,expected_baseline_modes)))


parse(p, varargin{:});
S = struct2cell(p.Results);
[baseline_fcn, baseline_mode, drug_abbr, ...
    epoch_length_seconds, epoch_overlap_seconds, p_threshold, show_epochs] = deal(S{:});

%Initialize frequency band names for labeling subplots
NC_band_names = {'\theta-waves (3-8 Hz)', '\alpha-waves (8-15 Hz)', '\beta-waves (15-30 Hz)'};
wPLI_band_names = {'\delta-waves (0.4-3 Hz)', '\theta-waves (3-8 Hz)', '\alpha-waves (8-15 Hz)', '\beta-waves (15-30 Hz)'};
PSD_band_names = {'\delta-waves (1-3 Hz)', '\theta-waves (3-8 Hz)', '\alpha-waves (8-15 Hz)', '\beta-waves (15-30 Hz)'};

if strcmp(feature, 'ASI')
    plot_counter = 0;
    h = figure;
    chs =  {[1, 2]
            [1, 3]
            [2, 4]
            [3, 4]
            [5, 6]};
    for ch_pair_i = 1:numel(chs)
        ch_pair = chs{ch_pair_i};
        [pre, post] = deal(cell(numel(predrug_features), 1));
        for sub = 1:numel(predrug_features)
                pre{sub} = predrug_features(sub).(feature)(ch_pair(1), ch_pair(2), :);
                post{sub} = postdrug_features(sub).(feature)(ch_pair(1), ch_pair(2), :);
        end
        ch_name = sprintf("%s vs %s", expected_channels{ch_pair(1)}, expected_channels{ch_pair(2)});
        plot_counter = plot_counter + 1;
        
        subplot(3, numel(chs), plot_counter);
        plot_timeline(pre, post, epoch_length_seconds, epoch_overlap_seconds, show_epochs);
        title(sprintf("%s, %s", feature, ch_name))
        if ch_pair_i == 1
            ylabel(feature, 'FontWeight','bold')
        end
        
        subplot(3, numel(chs), plot_counter+numel(chs));
        plot_bl_timeline(pre, post, ...
            epoch_length_seconds, epoch_overlap_seconds, ...
            show_epochs, ...
            baseline_fcn, baseline_mode);
        if strcmp(baseline_mode, 'relative')
            title(sprintf("%s/BL", drug_abbr))
        else
            title(sprintf("%s-BL", drug_abbr))
        end
        if ch_pair_i == 1
            ylabel(feature, 'FontWeight','bold')
        end
        xlabel('seconds')
        
        subplot(3, numel(chs), plot_counter+(2*numel(chs)));
        plot_bl_corr(pre, post, ...
            epoch_length_seconds, epoch_overlap_seconds, ...
            show_epochs, baseline_fcn, p_threshold);
        title(sprintf("Correlation %s-BL", ...
            drug_abbr));
        if ch_pair_i == 1
            ylabel("Corr coeff", 'FontWeight','bold')
        elseif ch_pair_i == numel(chs)
            yyaxis right;
            ylabel('p-value', 'FontWeight','bold');
            yyaxis left;
        end
        xlabel('seconds');
    end
    savefig(h, fullfile(results_path, 'ASI_all_plots.fig'));
    saveas(h, fullfile(results_path, 'ASI_all_plots.png'));
    
elseif strcmp(feature, 'SC')
    [pre, post] = deal(cell(numel(predrug_features), 1));
    for sub = 1:numel(predrug_features)
                pre{sub} = predrug_features(sub).(feature)';
                post{sub} = postdrug_features(sub).(feature)';
    end
    h = figure('Visible', 'Off');
    subplot(3, 1, 1);
    plot_timeline(pre, post, epoch_length_seconds, epoch_overlap_seconds, show_epochs);
    title(sprintf("%s timeline", feature))
    ylabel('Suppression')

    subplot(3, 1, 2);
    plot_bl_timeline(pre, post, ...
        epoch_length_seconds, epoch_overlap_seconds, ...
        show_epochs, ...
        baseline_fcn, baseline_mode);
    if strcmp(baseline_mode, 'relative')
        ylabel(sprintf("%s %s/BL", feature, drug_abbr))
    else
        ylabel(sprintf("%s %s-BL", feature, drug_abbr))
    end
    xlabel('seconds')
    
    subplot(3, 1, 3);
    plot_bl_corr(pre, post, ...
        epoch_length_seconds, epoch_overlap_seconds, ...
        show_epochs, ...
        baseline_fcn, p_threshold)
    title(sprintf("%s corr %s-BL", ...
                    feature, drug_abbr))

    xlabel('seconds');
    ylabel("Corr coeff", 'FontWeight','bold');
    yyaxis right;
    ylabel('p-value', 'FontWeight','bold');
    yyaxis left;
    savefig(h, fullfile(results_path, 'SC_all_plots.fig'));
    saveas(h, fullfile(results_path, 'SC_all_plots.png'));
    
elseif strcmp(feature, 'cPSD')
    plot_counter = 0;
    ax_timeline = figure('Visible', 'Off');
    ax_bl = figure('Visible', 'Off');
    ax_corr = figure('Visible', 'Off');
    bands = numel(PSD_band_names);
    chs =  {[1, 2]
            [1, 3]
            [2, 4]
            [3, 4]
            [5, 6]};
    for band = 1:bands
        for ch_pair_i = 1:5
            ch_pair = chs{ch_pair_i};
            [pre, post] = deal(cell(numel(predrug_features), 1));
            for sub = 1:numel(predrug_features)
                    pre{sub} = predrug_features(sub).(feature)(ch_pair(1), ch_pair(2), band, :);
                    post{sub} = postdrug_features(sub).(feature)(ch_pair(1), ch_pair(2), band, :);
            end

            ch_name = sprintf("%s vs %s", expected_channels{ch_pair(1)}, expected_channels{ch_pair(2)});
            band_name = PSD_band_names{band};
            plot_counter = plot_counter + 1;

            figure(ax_timeline)
            subplot(bands, size(chs, 1), plot_counter);
            plot_timeline(pre, post, epoch_length_seconds, epoch_overlap_seconds, show_epochs);
            if band == 1
                title(sprintf("%s, %s", feature, ch_name))
            end
            if ch_pair_i == 1
                ylabel(sprintf("%s", band_name), 'FontWeight','bold')
            end

            figure(ax_bl)
            subplot(bands, size(chs, 1), plot_counter);
            plot_bl_timeline(pre, post, ...
                epoch_length_seconds, epoch_overlap_seconds, ...
                show_epochs, ...
                baseline_fcn, baseline_mode);
            if band == 1
                if strcmp(baseline_mode, 'relative')
                    title(sprintf("%s %s/BL, %s", feature, drug_abbr, ch_name))
                else
                    title(sprintf("%s %s-BL, %s", feature, drug_abbr, ch_name))
                end
            end
            if ch_pair_i == 1
                ylabel(sprintf("%s", band_name), 'FontWeight','bold')
            end
            xlabel('seconds')

            figure(ax_corr)
            subplot(bands, size(chs, 1), plot_counter);
            plot_bl_corr(pre, post, ...
                epoch_length_seconds, epoch_overlap_seconds, ...
                show_epochs, baseline_fcn, p_threshold);
            if band == 1
                title(sprintf("%s %s corr %s-BL", ...
                    feature, ch_name, drug_abbr));
            elseif band == bands
                  xlabel('seconds');
            end
            if ch_pair_i == 1
                ylabel(sprintf("Corr coeff\n%s", band_name), 'FontWeight','bold')
            elseif ch_pair_i == 5
                yyaxis right;
                ylabel('p-value', 'FontWeight','bold')
                yyaxis left;
            end
        end
    end
    savefig(ax_timeline, fullfile(results_path, 'cPSD_timelines.fig'));
    savefig(ax_bl, fullfile(results_path, 'cPSD_baseline.fig'));
    savefig(ax_corr, fullfile(results_path, 'cPSD_correlation.fig'));
    saveas(ax_timeline, fullfile(results_path, 'cPSD_timelines.png'));
    saveas(ax_bl, fullfile(results_path, 'cPSD_baseline.png'));
    saveas(ax_corr, fullfile(results_path, 'cPSD_correlation.png'));
    
elseif strcmp(feature, 'wPLI')
    plot_counter = 0;
    ax_timeline = figure('Visible', 'Off');
    ax_bl = figure('Visible', 'Off');
    ax_corr = figure('Visible', 'Off');
    bands = numel(wPLI_band_names);
    chs =  {[1, 2]
            [1, 3]
            [2, 4]
            [3, 4]};
    for band = 1:bands
        for ch_pair_i = 1:numel(chs)
            ch_pair = chs{ch_pair_i};
            [pre, post] = deal(cell(numel(predrug_features), 1));
            for sub = 1:numel(predrug_features)
                    pre{sub} = predrug_features(sub).(feature)(ch_pair(1), ch_pair(2), band, :);
                    post{sub} = postdrug_features(sub).(feature)(ch_pair(1), ch_pair(2), band, :);
            end
            
            ch_name = sprintf("%s vs %s", expected_channels{ch_pair(1)}, expected_channels{ch_pair(2)});
            band_name = wPLI_band_names{band};
            plot_counter = plot_counter + 1;
            
            figure(ax_timeline)
            subplot(size(chs, 1), bands, plot_counter);
            plot_timeline(pre, post, epoch_length_seconds, epoch_overlap_seconds, show_epochs);
            if band == 1
                title(sprintf("%s, %s", feature, ch_name))
            elseif band == bands
                xlabel('seconds')
            end
            if ch_pair_i == 1
                ylabel(sprintf("%s", band_name), 'FontWeight','bold')
            end
            
            figure(ax_bl)
            subplot(size(chs, 1), bands, plot_counter);
            plot_bl_timeline(pre, post, ...
                epoch_length_seconds, epoch_overlap_seconds, ...
                show_epochs, ...
                baseline_fcn, baseline_mode);
            if band == 1
                if strcmp(baseline_mode, 'relative')
                    title(sprintf("%s %s/BL, %s", feature, drug_abbr, ch_name))
                else
                    title(sprintf("%s %s-BL, %s", feature, drug_abbr, ch_name))
                end
            elseif band == bands
                xlabel('seconds')
            end
            if ch_pair_i == 1
                ylabel(sprintf("%s", band_name), 'FontWeight','bold')
            end
            
            figure(ax_corr)
            subplot(size(chs, 1), bands, plot_counter);
            plot_bl_corr(pre, post, ...
                epoch_length_seconds, epoch_overlap_seconds, ...
                show_epochs, baseline_fcn, p_threshold);
            if band == 1
                title(sprintf("%s %s corr, %s-BL", ...
                    feature, ch_name, drug_abbr))
            elseif band == bands
                xlabel('seconds', 'FontWeight','bold')
            end 
            if ch_pair_i == 1
                ylabel(sprintf("Corr coeff\n%s", band_name), 'FontWeight','bold')
            elseif ch_pair_i == numel(chs)
                yyaxis right;
                ylabel('p-value', 'FontWeight','bold')
                yyaxis left;
            end
        end
    end   
    savefig(ax_timeline, fullfile(results_path, 'wPLI_timelines.fig'));
    savefig(ax_bl, fullfile(results_path, 'wPLI_baseline.fig'));
    savefig(ax_corr, fullfile(results_path, 'wPLI_correlation.fig'));
    saveas(ax_timeline, fullfile(results_path, 'wPLI_timelines.png'));
    saveas(ax_bl, fullfile(results_path, 'wPLI_baseline.png'));
    saveas(ax_corr, fullfile(results_path, 'wPLI_correlation.png'));
    
elseif strcmp(feature, 'NC') || strcmp(feature, 'PSD')
    plot_counter = 0;
    ax_timeline = figure('Visible', 'Off');
    ax_bl = figure('Visible', 'Off');
    ax_corr = figure('Visible', 'Off');
    if strcmp(feature, 'NC')
        band_names = NC_band_names;
        chs = 4;
    elseif strcmp(feature, 'PSD')
        band_names = PSD_band_names;
        chs = 8;
    end
    bands = numel(band_names);
    for band = 1:bands
        for ch = 1:chs
            
            [pre, post] = deal(cell(numel(predrug_features), 1));
            for sub = 1:numel(predrug_features)
                pre{sub} = predrug_features(sub).(feature)(ch, band, :);
                post{sub} = postdrug_features(sub).(feature)(ch, band, :);
            end
            ch_name = expected_channels{ch};
            plot_counter = plot_counter + 1;
            
            figure(ax_timeline)
            subplot(bands, chs, plot_counter);
            plot_timeline(pre, post, epoch_length_seconds, epoch_overlap_seconds, show_epochs);
            if band == 1
                title(sprintf("%s, %s", feature, ch_name))
            elseif band == bands
                xlabel('seconds')
            end
            if ch == 1
                ylabel(sprintf("%s", band_names{band}), 'FontWeight','bold')
            end
            
            figure(ax_bl)
            subplot(bands, chs, plot_counter);
            plot_bl_timeline(pre, post, ...
                epoch_length_seconds, epoch_overlap_seconds, ...
                show_epochs, ...
                baseline_fcn, baseline_mode);
            if band == 1
                if strcmp(baseline_mode, 'relative')
                    title(sprintf("%s %s/BL, %s", feature, drug_abbr, ch_name))
                else
                    title(sprintf("%s %s-BL, %s", feature, drug_abbr, ch_name))
                end
            elseif band == bands
                xlabel('seconds')
            end
            if ch == 1
                ylabel(sprintf("%s", band_names{band}), 'FontWeight','bold')
            end
            
            figure(ax_corr)
            subplot(bands, chs, plot_counter);
            plot_bl_corr(pre, post, ...
                epoch_length_seconds, epoch_overlap_seconds, ...
                show_epochs, baseline_fcn, p_threshold);
            if band == 1
                title(sprintf("%s %s corr, %s-BL", ...
                    feature, ch_name, drug_abbr), 'FontWeight','bold')
            elseif band == bands
                xlabel('seconds')
            end
            if ch == 1
                ylabel(sprintf("Corr coeff\n%s", band_names{band}), 'FontWeight','bold')
            elseif ch == chs
                yyaxis right;
                ylabel('p-value', 'FontWeight','bold')
                yyaxis left;
            end
        end
    end
    savefig(ax_timeline, fullfile(results_path, strcat(feature, '_timelines.fig')));
    savefig(ax_bl, fullfile(results_path, strcat(feature, '_baseline.fig')));
    savefig(ax_corr, fullfile(results_path, strcat(feature, '_correlation.fig')));
    saveas(ax_timeline, fullfile(results_path, strcat(feature, '_timelines.png')));
    saveas(ax_bl, fullfile(results_path, strcat(feature, '_baseline.png')));
    saveas(ax_corr, fullfile(results_path, strcat(feature, '_correlation.png')));
    
elseif any(strcmp(feature, {'aEEG', 'rEEG', 'MFDFA'}))
    subfeatures = fieldnames(predrug_features(1).(feature));
    for subfeature_i = 1:numel(subfeatures)
        h = figure('Visible', 'Off');
        plot_counter = 0;
        subfeature = subfeatures{subfeature_i};
        chs = 8;
        for ch = 1:chs
            [pre, post] = deal(cell(numel(predrug_features), 1));
            for sub = 1:numel(predrug_features)
                pre{sub} = predrug_features(sub).(feature).(subfeature)(ch, :)';
                post{sub} = postdrug_features(sub).(feature).(subfeature)(ch, :)';
            end
            ch_name = expected_channels{ch};
            plot_counter = plot_counter + 1;
            
            subplot(3, chs, plot_counter);
            plot_timeline(pre, post, epoch_length_seconds, epoch_overlap_seconds, show_epochs);
            title(sprintf("%s %s %s", feature, subfeature, ch_name))
            if ch == 1
                ylabel(sprintf("%s", subfeature), 'FontWeight','bold')
            end
            xlabel('seconds')

            subplot(3, chs, plot_counter+chs);
            plot_bl_timeline(pre, post, ...
                epoch_length_seconds, epoch_overlap_seconds, ...
                show_epochs, ...
                baseline_fcn, baseline_mode);
            if ch == 1
                if strcmp(baseline_mode, 'relative')
                    ylabel(sprintf("%s %s/BL", subfeature, drug_abbr), 'FontWeight','bold')
                else
                    ylabel(sprintf("%s %s-BL", subfeature, drug_abbr), 'FontWeight','bold')
                end
            end
            xlabel('seconds')
            
            subplot(3, chs, plot_counter+(chs*2));
            plot_bl_corr(pre, post, ...
                epoch_length_seconds, epoch_overlap_seconds, ...
                show_epochs, baseline_fcn, p_threshold);
            title(sprintf("Correlation, %s-BL", drug_abbr))
            if ch == 1
                ylabel(sprintf("%s Corr coeff", subfeature), 'FontWeight','bold')
            elseif ch == chs
                yyaxis right;
                ylabel('p-value', 'FontWeight','bold')
                yyaxis left;
            end
            xlabel('seconds')
        end
        savefig(h, fullfile(results_path, strcat(feature, '_', subfeature, '_all_plots.fig')));
        saveas(h, fullfile(results_path, strcat(feature, '_', subfeature, '_all_plots.png')))
    end
end
close all;
end