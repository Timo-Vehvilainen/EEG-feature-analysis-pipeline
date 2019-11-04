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
%           feature_names variable)
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
%           administration. Defaults to 5 predrug epochs and 10 postdrug epochs.
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
%           parameter determines whether the baseline timeline figures are
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
%       hide_p_values:
%           a boolean variable, indicating whether the correlation plots
%           want to be plotted as a double sided plot with correlation
%           coefficients on one side and p-values on the other (false), or
%           with just the correlation coefficient values (true).
%       corr_type:
%           a string variable indicating the type of correlation that is to
%           be used for calculating the coefficients. Has to have a value
%           of either 'Spearman', 'Pearson' or 'Kendall'

%Parse input:

feature_names = {...
    'ASI', ...
    'wPLI', ...
    'NC', ...
    'PSD', ...
    'cPSD', ...
    'aEEG', ...
    'rEEG', ...
    'MFDFA', ...
    'SC'};
correlation_types = {{'Spearman'}, {'Pearson'}, {'Kendall'}};

%check that given feature is valid.
if ~any(strcmp(feature_names, feature))
    error(sprintf("Unexpected feature: '%s'. Expected one of:\n", feature) + ...
        sprintf("%s ", feature_names{:}));
end

p = inputParser;
addParameter(p, 'epoch_length_seconds', 120, @isnumeric);
addParameter(p, 'epoch_overlap_seconds', 60, @isnumeric);
addParameter(p, 'show_epochs', [5, 10], @isnumeric);
addParameter(p, 'drug_abbr', 'DRUG', @ischar);
addParameter(p, 'p_threshold', 0.05, @isnumeric);
addParameter(p, 'hide_p_values', false, @islogical);

fcn_validation = @(x)  isa(x, 'function_handle');
addParameter(p, 'baseline_fcn', @baseline_median5, fcn_validation);

expected_baseline_modes = {'relative', 'absolute'};
addParameter(p, 'baseline_mode', 'relative', ...
    @(x) any(validatestring(x,expected_baseline_modes)));

addParameter(p, 'corr_type', 'Spearman', ...
    @(x) any(validatestring(x, [correlation_types{:}])));

parse(p, varargin{:});
S = struct2cell(p.Results);
[baseline_fcn, baseline_mode, corr_type, drug_abbr, ...
    epoch_length_seconds, epoch_overlap_seconds, ...
    hide_p_values, p_threshold, show_epochs] = deal(S{:});

wilcoxon_baseline_fcn = @baseline_last_epoch;

%Channel label names
channel_names = {'F3', 'F4', 'P3', 'P4', ...
    'Left', 'Right', 'Frontal', 'Parietal'};

%Initialize frequency band names for labeling subplots
NC_band_names = {'\theta-waves (3-8 Hz)', '\alpha-waves (8-15 Hz)', ...
    '\beta-waves (15-30 Hz)'};
wPLI_band_names = {'\delta-waves (0.4-3 Hz)', '\theta-waves (3-8 Hz)', ...
    '\alpha-waves (8-15 Hz)', '\beta-waves (15-30 Hz)'};
PSD_band_names = {'\delta-waves (1-3 Hz)', '\theta-waves (3-8 Hz)', ...
    '\alpha-waves (8-15 Hz)', '\beta-waves (15-30 Hz)'};

%Initialize correlation coefficient label based on the chose type of correlation
corr_labels = ["Spearman's \rho", "Pearson's r", "Kendall's \tau"];
corr_label_idx = strcmp([correlation_types{:}], corr_type);
corr_label = corr_labels(corr_label_idx);

%Choose a code block to handle plotting depenging on the feature:
% (the code blocks are similar in structure but different enough, such that 
% they have been chosen to be written separately instead of trying to
% generalize the code for all of them. If you want to add a new feature to
% the list of features to me plotted, you need to write a new block in this
% if-elseif decision tree and add the feature to the feature_names variable)
if strcmp(feature, 'ASI')
    plot_counter = 0;   %variable for counting subplots
    h = figure('Visible', 'Off', 'Position', get(0, 'Screensize'));
    
    %List of channels pairs to be plotted against each other. The indices
    %refer to the channels in the variable channel_names:
    chs =  {[1, 2]      %F3 vs F4
            [1, 3]      %F3 vs P3
            [2, 4]      %F4 vs P4
            [3, 4]      %P3 vs P4
            [5, 6]};    %Left vs Right
        
    %initialize the data matrix that is to be saved with the plots
    data_mat = cell(4, numel(chs));
    
    %go through the channel pairs
    for ch_pair_i = 1:numel(chs)
        ch_pair = chs{ch_pair_i};
        
        %Extract to predrug and postdrug features for the selected channels
        [pre, post] = deal(cell(numel(predrug_features), 1));
        for sub = 1:numel(predrug_features)
                pre{sub} = predrug_features(sub).(feature)(ch_pair(1), ch_pair(2), :);
                post{sub} = postdrug_features(sub).(feature)(ch_pair(1), ch_pair(2), :);
        end
        
        %initialize channel name for titles
        ch_name = sprintf("%s vs %s", channel_names{ch_pair(1)}, channel_names{ch_pair(2)});
        plot_counter = plot_counter + 1;
        
        %Plot absolute timeline
        subplot(4, numel(chs), plot_counter);
        data_mat{1, ch_pair_i} = plot_timeline(pre, post, epoch_length_seconds, epoch_overlap_seconds, show_epochs);
        %Title & labels
        title(sprintf("%s, %s", feature, ch_name))
        if ch_pair_i == 1
            ylabel(feature, 'FontWeight','bold')
        end
        
        %Plot signrank
        subplot(4, numel(chs), plot_counter+numel(chs));
        data_mat{2, ch_pair_i} = plot_signrank(pre, post, ...
            epoch_length_seconds, epoch_overlap_seconds, ...
            show_epochs, wilcoxon_baseline_fcn, p_threshold, hide_p_values);
        %Title & labels
        title(sprintf("Signrank %s-BL", ...
            drug_abbr));
        if ch_pair_i == 1
            ylabel(sprintf("median(%s-BL)", drug_abbr), 'FontWeight','bold')
        elseif ch_pair_i == numel(chs) && ~hide_p_values
            yyaxis right;
            ylabel('p-value', 'FontWeight','bold');
            yyaxis left;
        end
        xlabel('seconds');
        
        %Plot baseline-relative timeline
        subplot(4, numel(chs), plot_counter+(2*numel(chs)));
        data_mat{3, ch_pair_i} = plot_bl_timeline(pre, post, ...
            epoch_length_seconds, epoch_overlap_seconds, ...
            show_epochs, ...
            baseline_fcn, baseline_mode);
        %Title & labels
        if strcmp(baseline_mode, 'relative')
            title(sprintf("%s/BL", drug_abbr))
        else
            title(sprintf("%s-BL", drug_abbr))
        end
        if ch_pair_i == 1
            ylabel(feature, 'FontWeight','bold')
        end
        xlabel('seconds')
        
        %Plot baseline-delta correlation
        subplot(4, numel(chs), plot_counter+(3*numel(chs)));
        data_mat{4, ch_pair_i} = plot_bl_delta_corr(pre, post, ...
            epoch_length_seconds, epoch_overlap_seconds, ...
            show_epochs, baseline_fcn, p_threshold, 'delta', hide_p_values, corr_type);
        %Title & labels
        title(sprintf("Correlation \\Delta%s-BL", ...
            drug_abbr));
        if ch_pair_i == 1
            ylabel(corr_label, 'FontWeight','bold')
        elseif ch_pair_i == numel(chs) && ~hide_p_values
            yyaxis right;
            ylabel('p-value', 'FontWeight','bold');
            yyaxis left;
        end
        xlabel('seconds');
    end
    %Save the figure and data matrix in the results-path
    savefig(h, fullfile(results_path, 'ASI_all_plots.fig'));
    export_fig(h, fullfile(results_path, 'ASI_all_plots.pdf'), '-transparent');
    save(fullfile(results_path, 'ASI_all_plots.mat'), 'data_mat');
    
elseif strcmp(feature, 'SC')
    data_mat = cell(4, 1);
    
    %Extract to predrug and postdrug feature for SC
    [pre, post] = deal(cell(numel(predrug_features), 1));
    for sub = 1:numel(predrug_features)
                pre{sub} = predrug_features(sub).(feature)';
                post{sub} = postdrug_features(sub).(feature)';
    end
    h = figure('Visible', 'Off');
    
    %Plot absolute timeline
    subplot(4, 1, 1);
    data_mat{1, 1} = plot_timeline(pre, post, ...
            epoch_length_seconds, epoch_overlap_seconds, show_epochs);
    title(sprintf("%s timeline", feature))
    ylabel('Suppression')
    
    %Plot signrank
    subplot(4, 1, 2);
    data_mat{2, 1} = plot_signrank(pre, post, ...
            epoch_length_seconds, epoch_overlap_seconds, ...
            show_epochs, wilcoxon_baseline_fcn, p_threshold, hide_p_values);
    %Title & labels
    title(sprintf("Signrank %s-BL", ...
        drug_abbr));
    ylabel(sprintf("median(%s-BL)", drug_abbr), 'FontWeight','bold');
    if ~hide_p_values
        yyaxis right;
        ylabel('p-value', 'FontWeight','bold');
        yyaxis left;
    end
    xlabel('seconds');

    %Plot baseline-relative timeline
    subplot(4, 1, 3);
    data_mat{3, 1} = plot_bl_timeline(pre, post, ...
        epoch_length_seconds, epoch_overlap_seconds, ...
        show_epochs, ...
        baseline_fcn, baseline_mode);
    %Title & labels
    if strcmp(baseline_mode, 'relative')
        title(sprintf("%s/BL", drug_abbr))
    else
        title(sprintf("%s-BL", drug_abbr))
    end
    if strcmp(baseline_mode, 'relative')
        ylabel(sprintf("%s %s/BL", feature, drug_abbr))
    else
        ylabel(sprintf("%s %s-BL", feature, drug_abbr))
    end
    xlabel('seconds')
    
    %Plot baseline-delta correlation
    subplot(4, 1, 4);
    data_mat{4, 1} = plot_bl_delta_corr(pre, post, ...
        epoch_length_seconds, epoch_overlap_seconds, ...
        show_epochs, ...
        baseline_fcn, p_threshold, 'delta', hide_p_values, corr_type);
    %Title & labels
    title(sprintf("%s corr \\Delta%s-BL", ...
                    feature, drug_abbr))
    xlabel('seconds');
    ylabel(corr_label, 'FontWeight','bold');
    if ~hide_p_values
        yyaxis right;
        ylabel('p-value', 'FontWeight','bold');
        yyaxis left;
    end
    
    savefig(h, fullfile(results_path, 'SC_all_plots.fig'));
    export_fig(h, fullfile(results_path, 'SC_all_plots.pdf'), '-transparent');
    save(fullfile(results_path, 'SC_all_plots.mat'), 'data_mat');
    
elseif strcmp(feature, 'cPSD')
    plot_counter = 0;
    %Plot the statistics in separate figures
    ax_timeline = figure('Visible', 'Off', 'Position', get(0, 'Screensize'));
    ax_signrank = figure('Visible', 'Off', 'Position', get(0, 'Screensize'));
    ax_bl = figure('Visible', 'Off', 'Position', get(0, 'Screensize'));
    ax_corr = figure('Visible', 'Off', 'Position', get(0, 'Screensize'));
    %Number of frequency bands
    bands = numel(PSD_band_names);
    %Channel pairs
    chs =  {[1, 2]
            [1, 3]
            [2, 4]
            [3, 4]
            [5, 6]};
        
    %Plotted data matrices to be exported
    timeline_data_mat = cell(bands, numel(chs));
    signrank_data_mat = cell(bands, numel(chs));
    bl_data_mat = cell(bands, numel(chs));
    corr_data_mat = cell(bands, numel(chs));
    for band = 1:bands
        for ch_pair_i = 1:5
            ch_pair = chs{ch_pair_i};
            
            %Extract to predrug and postdrug features for the selected channels
            [pre, post] = deal(cell(numel(predrug_features), 1));
            for sub = 1:numel(predrug_features)
                    pre{sub} = predrug_features(sub).(feature)(ch_pair(1), ch_pair(2), band, :);
                    post{sub} = postdrug_features(sub).(feature)(ch_pair(1), ch_pair(2), band, :);
            end

            ch_name = sprintf("%s vs %s", channel_names{ch_pair(1)}, channel_names{ch_pair(2)});
            band_name = PSD_band_names{band};
            plot_counter = plot_counter + 1;

            %Plot absolute timeline
            figure(ax_timeline)
            subplot(bands, size(chs, 1), plot_counter);
            timeline_data_mat{band, ch_pair_i} = plot_timeline(pre, post, epoch_length_seconds, epoch_overlap_seconds, show_epochs);
            %Title & labels
            if band == 1
                title(sprintf("%s, %s", feature, ch_name))
            end
            if ch_pair_i == 1
                ylabel(sprintf("%s", band_name), 'FontWeight','bold')
            end
            
            %Plot signrank
            figure(ax_signrank)
            subplot(bands, size(chs, 1), plot_counter);
            signrank_data_mat{band, ch_pair_i} = plot_signrank(pre, post, ...
                epoch_length_seconds, epoch_overlap_seconds, ...
                show_epochs, wilcoxon_baseline_fcn, p_threshold, hide_p_values);
            %Title & labels
            if band == 1
                title(sprintf("%s %s Signrank %s-BL", ...
                    feature, ch_name, drug_abbr));
            elseif band == bands
                  xlabel('seconds');
            end
            if ch_pair_i == 1
                ylabel(sprintf("median(%s-BL)\n%s", drug_abbr, band_name), 'FontWeight','bold')
            elseif ch_pair_i == 5 && ~hide_p_values
                yyaxis right;
                ylabel('p-value', 'FontWeight','bold')
                yyaxis left;
            end

            %Plot baseline-relative timeline
            figure(ax_bl)
            subplot(bands, size(chs, 1), plot_counter);
            bl_data_mat{band, ch_pair_i} = plot_bl_timeline(pre, post, ...
                epoch_length_seconds, epoch_overlap_seconds, ...
                show_epochs, ...
                baseline_fcn, baseline_mode);
            %Title & labels
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

            %Plot baseline-delta correlation
            figure(ax_corr)
            subplot(bands, size(chs, 1), plot_counter);
            corr_data_mat{band, ch_pair_i} = plot_bl_delta_corr(pre, post, ...
                epoch_length_seconds, epoch_overlap_seconds, ...
                show_epochs, baseline_fcn, p_threshold, 'delta', hide_p_values, corr_type);
            %Title & labels
            if band == 1
                title(sprintf("%s %s corr \\Delta%s-BL", ...
                    feature, ch_name, drug_abbr));
            elseif band == bands
                  xlabel('seconds');
            end
            if ch_pair_i == 1
                ylabel(sprintf("%s\n%s", corr_label, band_name), 'FontWeight','bold')
            elseif ch_pair_i == 5 && ~hide_p_values
                yyaxis right;
                ylabel('p-value', 'FontWeight','bold')
                yyaxis left;
            end
        end
    end
    %Save all figures and matrices
    savefig(ax_timeline, fullfile(results_path, 'cPSD_timelines.fig'));
    savefig(ax_signrank, fullfile(results_path, 'cPSD_signrank.fig'));
    savefig(ax_bl, fullfile(results_path, 'cPSD_baseline.fig'));
    savefig(ax_corr, fullfile(results_path, 'cPSD_correlation_delta.fig'));
    
    export_fig(ax_timeline, fullfile(results_path, 'cPSD_timelines.pdf'), '-transparent');
    export_fig(ax_signrank, fullfile(results_path, 'cPSD_signrank.pdf'), '-transparent');
    export_fig(ax_bl, fullfile(results_path, 'cPSD_baseline.pdf'), '-transparent');
    export_fig(ax_corr, fullfile(results_path, 'cPSD_correlation_delta.pdf'), '-transparent');
    
    save(fullfile(results_path, 'cPSD_timelines.mat'), 'timeline_data_mat');
    save(fullfile(results_path, 'cPSD_signrank.mat'), 'signrank_data_mat');
    save(fullfile(results_path, 'cPSD_baseline.mat'), 'bl_data_mat');
    save(fullfile(results_path, 'cPSD_correlation_delta.mat'), 'corr_data_mat');
    
elseif strcmp(feature, 'wPLI')
    plot_counter = 0;
    ax_timeline = figure('Visible', 'Off', 'Position', get(0, 'Screensize'));
    ax_signrank = figure('Visible', 'Off', 'Position', get(0, 'Screensize'));
    ax_bl = figure('Visible', 'Off', 'Position', get(0, 'Screensize'));
    ax_corr = figure('Visible', 'Off', 'Position', get(0, 'Screensize'));
    bands = numel(wPLI_band_names);
    chs =  {[1, 2]
            [1, 3]
            [2, 4]
            [3, 4]};
    timeline_data_mat = cell(bands, numel(chs));
    signrank_data_mat = cell(bands, numel(chs));
    bl_data_mat = cell(bands, numel(chs));
    corr_data_mat = cell(bands, numel(chs));
    for band = 1:bands
        for ch_pair_i = 1:numel(chs)
            ch_pair = chs{ch_pair_i};
            
            %Extract to predrug and postdrug features for the selected channels
            [pre, post] = deal(cell(numel(predrug_features), 1));
            for sub = 1:numel(predrug_features)
                    pre{sub} = predrug_features(sub).(feature)(ch_pair(1), ch_pair(2), band, :);
                    post{sub} = postdrug_features(sub).(feature)(ch_pair(1), ch_pair(2), band, :);
            end
            
            ch_name = sprintf("%s vs %s", channel_names{ch_pair(1)}, channel_names{ch_pair(2)});
            band_name = wPLI_band_names{band};
            plot_counter = plot_counter + 1;
            
            %Plot absolute timeline
            figure(ax_timeline)
            subplot(size(chs, 1), bands, plot_counter);
            timeline_data_mat{band, ch_pair_i} = plot_timeline(pre, post, epoch_length_seconds, epoch_overlap_seconds, show_epochs);
            %Title & labels
            if band == 1
                title(sprintf("%s, %s", feature, ch_name))
            elseif band == bands
                xlabel('seconds')
            end
            if ch_pair_i == 1
                ylabel(sprintf("%s", band_name), 'FontWeight','bold')
            end
            
            %Plot signrank
            figure(ax_signrank)
            subplot(size(chs, 1), bands, plot_counter);
            signrank_data_mat{band, ch_pair_i} = plot_signrank(pre, post, ...
                epoch_length_seconds, epoch_overlap_seconds, ...
                show_epochs, wilcoxon_baseline_fcn, p_threshold, hide_p_values);
            %Title & labels
            if band == 1
                title(sprintf("%s %s Signrank, %s-BL", ...
                    feature, ch_name, drug_abbr))
            elseif band == bands
                xlabel('seconds', 'FontWeight','bold')
            end 
            if ch_pair_i == 1
                ylabel(sprintf("median(%s-BL)\n%s", drug_abbr, band_name), 'FontWeight','bold')
            elseif ch_pair_i == numel(chs) && ~hide_p_values
                yyaxis right;
                ylabel('p-value', 'FontWeight','bold')
                yyaxis left;
            end
            
            %Plot baseline-relative timeline
            figure(ax_bl)
            subplot(size(chs, 1), bands, plot_counter);
            bl_data_mat{band, ch_pair_i} = plot_bl_timeline(pre, post, ...
                epoch_length_seconds, epoch_overlap_seconds, ...
                show_epochs, ...
                baseline_fcn, baseline_mode);
            %Title & labels
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
            
            %Plot baseline-delta correlation
            figure(ax_corr)
            subplot(size(chs, 1), bands, plot_counter);
            corr_data_mat{band, ch_pair_i} = plot_bl_delta_corr(pre, post, ...
                epoch_length_seconds, epoch_overlap_seconds, ...
                show_epochs, baseline_fcn, p_threshold, 'delta', hide_p_values, corr_type);
            %Title & labels
            if band == 1
                title(sprintf("%s %s corr, \\Delta%s-BL", ...
                    feature, ch_name, drug_abbr))
            elseif band == bands
                xlabel('seconds', 'FontWeight','bold')
            end 
            if ch_pair_i == 1
                ylabel(sprintf("%s\n%s", corr_label, band_name), 'FontWeight','bold')
            elseif ch_pair_i == numel(chs) && ~hide_p_values
                yyaxis right;
                ylabel('p-value', 'FontWeight','bold')
                yyaxis left;
            end
        end
    end   
    %Save all figures and matrices
    savefig(ax_timeline, fullfile(results_path, 'wPLI_timelines.fig'));
    savefig(ax_signrank, fullfile(results_path, 'wPLI_signrank.fig'));
    savefig(ax_bl, fullfile(results_path, 'wPLI_baseline.fig'));
    savefig(ax_corr, fullfile(results_path, 'wPLI_correlation_delta.fig'));
    
    export_fig(ax_timeline, fullfile(results_path, 'wPLI_timelines.pdf'), '-transparent');
    export_fig(ax_signrank, fullfile(results_path, 'wPLI_signrank.pdf'), '-transparent');
    export_fig(ax_bl, fullfile(results_path, 'wPLI_baseline.pdf'), '-transparent');
    export_fig(ax_corr, fullfile(results_path, 'wPLI_correlation_delta.pdf'), '-transparent');
    
    save(fullfile(results_path, 'wPLI_timelines.mat'), 'timeline_data_mat');
    save(fullfile(results_path, 'wPLI_signrank.mat'), 'signrank_data_mat');
    save(fullfile(results_path, 'wPLI_baseline.mat'), 'bl_data_mat');
    save(fullfile(results_path, 'wPLI_correlation_delta.mat'), 'corr_data_mat');
    
%Handling NC and PSD is similar enough to be combined into a single block:
elseif strcmp(feature, 'NC') || strcmp(feature, 'PSD')
    plot_counter = 0;
    ax_timeline = figure('Visible', 'Off', 'Position', get(0, 'Screensize'));
    ax_signrank = figure('Visible', 'Off', 'Position', get(0, 'Screensize'));
    ax_bl = figure('Visible', 'Off', 'Position', get(0, 'Screensize'));
    ax_corr = figure('Visible', 'Off', 'Position', get(0, 'Screensize'));
    
    %Choose frequency bands depending on the feature
    if strcmp(feature, 'NC')
        band_names = NC_band_names;
        chs = 4;    %number of channels to plot
    elseif strcmp(feature, 'PSD')
        band_names = PSD_band_names;
        chs = 8;
    end
    bands = numel(band_names);
    
    timeline_data_mat = cell(bands, chs);
    signrank_data_mat = cell(bands, chs);
    bl_data_mat = cell(bands, chs);
    corr_data_mat = cell(bands, chs);
    %loop through frequency bands
    for band = 1:bands
        %loop through channels
        for ch = 1:chs
            
            %Extract to predrug and postdrug features for the selected channels
            [pre, post] = deal(cell(numel(predrug_features), 1));
            for sub = 1:numel(predrug_features)
                pre{sub} = predrug_features(sub).(feature)(ch, band, :);
                post{sub} = postdrug_features(sub).(feature)(ch, band, :);
            end
            
            ch_name = channel_names{ch};
            plot_counter = plot_counter + 1;
            
            %Plot absolute timeline
            figure(ax_timeline)
            subplot(bands, chs, plot_counter);
            timeline_data_mat{band, ch} = plot_timeline(pre, post, epoch_length_seconds, epoch_overlap_seconds, show_epochs);
            %Title & labels
            if band == 1
                title(sprintf("%s, %s", feature, ch_name))
            elseif band == bands
                xlabel('seconds')
            end
            if ch == 1
                ylabel(sprintf("%s", band_names{band}), 'FontWeight','bold')
            end
            
            %Plot signrank
            figure(ax_signrank)
            subplot(bands, chs, plot_counter);
            signrank_data_mat{band, ch} = plot_signrank(pre, post, ...
                epoch_length_seconds, epoch_overlap_seconds, ...
                show_epochs, wilcoxon_baseline_fcn, p_threshold, hide_p_values);
            %Title & labels
            if band == 1
                title(sprintf("%s %s Signrank, %s-BL", ...
                    feature, ch_name, drug_abbr), 'FontWeight','bold')
            elseif band == bands
                xlabel('seconds')
            end
            if ch == 1
                ylabel(sprintf("median(%s-BL)\n%s", drug_abbr, band_names{band}), 'FontWeight','bold')
            elseif ch == chs && ~hide_p_values
                yyaxis right;
                ylabel('p-value', 'FontWeight','bold')
                yyaxis left;
            end
            
            %Plot baseline-relative timeline
            figure(ax_bl)
            subplot(bands, chs, plot_counter);
            bl_data_mat{band, ch} = plot_bl_timeline(pre, post, ...
                epoch_length_seconds, epoch_overlap_seconds, ...
                show_epochs, ...
                baseline_fcn, baseline_mode);
            %Title & labels
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
            
            %Plot baseline-delta correlation
            figure(ax_corr)
            subplot(bands, chs, plot_counter);
            corr_data_mat{band, ch} = plot_bl_delta_corr(pre, post, ...
                epoch_length_seconds, epoch_overlap_seconds, ...
                show_epochs, baseline_fcn, p_threshold, 'delta', hide_p_values, corr_type);
            %Title & labels
            if band == 1
                title(sprintf("%s %s corr, \\Delta%s-BL", ...
                    feature, ch_name, drug_abbr), 'FontWeight','bold')
            elseif band == bands
                xlabel('seconds')
            end
            if ch == 1
                ylabel(sprintf("%s\n%s", corr_label, band_names{band}), 'FontWeight','bold')
            elseif ch == chs && ~hide_p_values
                yyaxis right;
                ylabel('p-value', 'FontWeight','bold')
                yyaxis left;
            end
        end
    end
    %Save all figures and matrices
    savefig(ax_timeline, fullfile(results_path, strcat(feature, '_timelines.fig')));
    savefig(ax_signrank, fullfile(results_path, strcat(feature, '_signrank.fig')));
    savefig(ax_bl, fullfile(results_path, strcat(feature, '_baseline.fig')));
    savefig(ax_corr, fullfile(results_path, strcat(feature, '_correlation_delta.fig')));
    
    export_fig(ax_timeline, fullfile(results_path, strcat(feature, '_timelines.pdf')), '-transparent');
    export_fig(ax_signrank, fullfile(results_path, strcat(feature, '_signrank.pdf')), '-transparent');
    export_fig(ax_bl, fullfile(results_path, strcat(feature, '_baseline.pdf')), '-transparent');
    export_fig(ax_corr, fullfile(results_path, strcat(feature, '_correlation_delta.pdf')), '-transparent');
    
    save(fullfile(results_path, strcat(feature, '_timelines.mat')), 'timeline_data_mat');
    save(fullfile(results_path, strcat(feature, '_signrank.mat')), 'signrank_data_mat');
    save(fullfile(results_path, strcat(feature, '_baseline.mat')), 'bl_data_mat');
    save(fullfile(results_path, strcat(feature, '_correlation_delta.mat')), 'corr_data_mat');
%Handling aEEG, rEEG and MFDFA with one code block:
elseif any(strcmp(feature, {'aEEG', 'rEEG', 'MFDFA'}))
    
    %Handle each subfeature separately
    subfeatures = fieldnames(predrug_features(1).(feature));
    for subfeature_i = 1:numel(subfeatures)
        h = figure('Visible', 'Off', 'Position', get(0, 'Screensize'));
        plot_counter = 0;
        subfeature = subfeatures{subfeature_i};
        chs = 8;
        data_mat = cell(4, chs);
        
        for ch = 1:chs
            
            %Extract to predrug and postdrug features for the selected channels
            [pre, post] = deal(cell(numel(predrug_features), 1));
            for sub = 1:numel(predrug_features)
                pre{sub} = predrug_features(sub).(feature).(subfeature)(ch, :)';
                post{sub} = postdrug_features(sub).(feature).(subfeature)(ch, :)';
            end
            ch_name = channel_names{ch};
            plot_counter = plot_counter + 1;
            
            %Plot absolute timeline
            subplot(4, chs, plot_counter);
            data_mat{1, ch} = plot_timeline(pre, post, epoch_length_seconds, epoch_overlap_seconds, show_epochs);
            %Title & labels
            title(sprintf("%s %s %s", feature, subfeature, ch_name))
            if ch == 1
                ylabel(sprintf("%s", subfeature), 'FontWeight','bold')
            end
            xlabel('seconds')
            
            %Plot signrank
            subplot(4, chs, plot_counter+chs);
            data_mat{2, ch} = plot_signrank(pre, post, ...
                epoch_length_seconds, epoch_overlap_seconds, ...
                show_epochs, wilcoxon_baseline_fcn, p_threshold, hide_p_values);
            %Title & labels
            title(sprintf("Signrank, %s-BL", drug_abbr))
            if ch == 1
                ylabel(sprintf("%s median(%s-BL)", subfeature, drug_abbr), 'FontWeight','bold')
            elseif ch == chs && ~hide_p_values
                yyaxis right;
                ylabel('p-value', 'FontWeight','bold')
                yyaxis left;
            end
            xlabel('seconds')

            %Plot baseline-relative timeline
            subplot(4, chs, plot_counter+(chs*2));
            data_mat{3, ch} = plot_bl_timeline(pre, post, ...
                epoch_length_seconds, epoch_overlap_seconds, ...
                show_epochs, ...
                baseline_fcn, baseline_mode);
            %Title & labels
            if strcmp(baseline_mode, 'relative')
                title(sprintf("%s/BL", drug_abbr))
            else
                title(sprintf("%s-BL", drug_abbr))
            end
            if ch == 1
                if strcmp(baseline_mode, 'relative')
                    ylabel(sprintf("%s %s/BL", subfeature, drug_abbr), 'FontWeight','bold')
                else
                    ylabel(sprintf("%s %s-BL", subfeature, drug_abbr), 'FontWeight','bold')
                end
            end
            xlabel('seconds')
            
            %Plot baseline-delta correlation
            subplot(4, chs, plot_counter+(chs*3));
            data_mat{4, ch} = plot_bl_delta_corr(pre, post, ...
                epoch_length_seconds, epoch_overlap_seconds, ...
                show_epochs, baseline_fcn, p_threshold, 'delta', hide_p_values, corr_type);
            %Title & labels
            title(sprintf("Correlation, \\Delta%s-BL", drug_abbr))
            if ch == 1
                ylabel(sprintf("%s %s", subfeature, corr_label), 'FontWeight','bold')
            elseif ch == chs && ~hide_p_values
                yyaxis right;
                ylabel('p-value', 'FontWeight','bold')
                yyaxis left;
            end
            xlabel('seconds')
        end
        %Save all figures and matrices
        savefig(h, fullfile(results_path, strcat(feature, '_', subfeature, '_all_plots.fig')));
        export_fig(h, fullfile(results_path, strcat(feature, '_', subfeature, '_all_plots.pdf')), '-transparent');
        save(fullfile(results_path, strcat(feature, '_', subfeature, '_all_plots.mat')), 'data_mat');
    end
end
close all;
end