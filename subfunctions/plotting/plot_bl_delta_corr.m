%Timo Vehviläinen August 2019

function [data_table] = plot_bl_delta_corr(pre, post, epoch_length_seconds, epoch_overlap_seconds, ...
    show_epochs, baseline_fcn, p_threshold, correlation_mode, hide_p_values, correlation_type)
%PLOT_BL_CORR 
% This function takes in epochs of pre- & postdrug EEG data for a single
% feature, calculates a baseline from the predrug epochs, calculates
% correlation coeffiecients for each postdrug epoch against said baseline
% (and their p-values), and finally plots them in a two-sided graph.
%PARAMETERS:
%       pre:
%           A cell array, where each cell contains a vector or numeric
%           data, representing predrug feature values of epoched EEG data.
%       post:
%           A cell array, where each cell contains a vector or numeric
%           data, representing postdrug feature values of epoched EEG data.
%       epoch_length_seconds:
%           the duration of one epoch in seconds. This is used fro scaling
%           the x-axis.
%       epoch_overlap_seconds:
%           The amount of seconds that consecutive epochs overlap with one
%           another. This is also used for calculating the proper scaling
%           for the x-axis of the graph.
%       show_epochs:
%           A two-element numeric vector containing positive integers,
%           representing how many predrug & postdrug epochs (respectively)
%           are to be plotted from the data. (only the second value is
%           used.)
%       baseline_fcn:
%           A function variable, which takes a single numeric vector as
%           it's argument. This funciton is used for calculating the
%           baseline from a single predrug feature vector.
%       p_threshold:
%           a value between 0 and 1, representing the significance
%           threshold chosen for the correlations. Correlations where the
%           p-value is smaller than p_threshold will be highlighted in the
%           graph.
%       correlation_mode:
%           a string variable, used for determining whether the
%           correlations are calculated from [BL vs (DRUG-BL)] ('delta') or from 
%           [BL vs DRUG] (any other string)
%       hide_p_values:
%           a boolean variable, indicating whether the correlation plots
%           want to be plotted as a double sided plot with correlation
%           coefficients on one side and p-values on the other (false), or
%           with just the correlation coefficient values (true).
%       corr_type:
%           a string variable indicating the type of correlation that is to
%           be used for calculating the coefficients. Has to have a value
%           of either 'Spearman', 'Pearson' or 'Kendall'
%OUTPUT
%   Returns a table of column variables, where the first column variable
%   holds the x-axis values for the plot. The second variable holds the
%   correlation coefficients for each epoch, and the third variable holds
%   their p-values.

max_postdrug_epochs = min(show_epochs(2), max(cellfun(@numel, post)));
post_x = (epoch_length_seconds/2): ...
    (epoch_length_seconds - epoch_overlap_seconds): ...
    ((epoch_length_seconds*(max_postdrug_epochs)/2));

data_table = table(post_x', 'VariableNames', "X_axis");
all_baselines = cellfun(@(x)baseline_fcn(x), pre);
[ccs, ps] = deal(zeros(1, max_postdrug_epochs));

for epoch = 1:max_postdrug_epochs

    if strcmp(correlation_mode, 'delta')
        post_values = cellfun(@(x)x(epoch), post) - all_baselines;
    else
        post_values = cellfun(@(x)x(epoch), post);
    end

    [CC, P] = corr([post_values, all_baselines], ...
        'rows','pairwise', 'Type', correlation_type);
    ccs(epoch) = CC(1, 2);
    ps(epoch) = P(1, 2);
end

data_table = addvars(data_table, ccs', ps', 'NewVariableNames', {'CorrCoeffs', 'P_values'});

p_significant = (ps < p_threshold);

%Plot a horizontal line at the chosen significance level
epochs = 1:max_postdrug_epochs;

% separate the correlation coefficients and p-vals into 
% significant and unsignificant sets
sign_p_points = ps;
sign_p_points(~p_significant) = NaN;
unsign_p_points = ps;
unsign_p_points(p_significant) = NaN;

sign_cc_points = ccs;
sign_cc_points(~p_significant) = NaN;
unsign_cc_points = ccs;
unsign_cc_points(p_significant) = NaN;

if ~hide_p_values
    yyaxis right
    plot([0, max(post_x)], [p_threshold, p_threshold], '--');
    hold on;

    %Plot the p-values (differently depending on if they are below the
    %significance level)
    h1 = plot(post_x, unsign_p_points, 'rx', 'LineWidth', 1, 'MarkerSize', 10);
    hold on
    h2 = plot(post_x, sign_p_points, 'ro', 'LineWidth', 3, 'MarkerSize', 10);
    ylim([0, 1])
    hold off;

    %insert data tips
    enableDefaultInteractivity(gca);
    if ~isempty(h1)
        dtt = h1.DataTipTemplate;
        dtt.DataTipRows(1) = dataTipTextRow('p =', unsign_p_points);
        dtt.DataTipRows(2) = dataTipTextRow('postdrug epoch', epochs);
    end
    if ~isempty(h2)
        dtt = h2.DataTipTemplate;
        dtt.DataTipRows(1) = dataTipTextRow('p =', sign_p_points);
        dtt.DataTipRows(2) = dataTipTextRow('postdrug epoch', epochs);
    end

    %Plot the correlation coefficients (differently depending on if they have a
    %significant p-value or not).
    yyaxis left

    h3 = plot(post_x, unsign_cc_points, 'bx', 'MarkerSize', 10, 'LineWidth', 1);
    hold on;
    h4 = plot(post_x, sign_cc_points, 'bo', 'MarkerSize', 10, 'LineWidth', 3);
    hold off;
    ylim([-1, 1])

    %Insert data tips
    enableDefaultInteractivity(gca);
    if ~isempty(h3)
        dtt = h3.DataTipTemplate;
        dtt.DataTipRows(1) = dataTipTextRow('r = ', unsign_cc_points);
        dtt.DataTipRows(2) = dataTipTextRow('postdrug epoch', epochs);
    end
    if ~isempty(h4)
        dtt = h4.DataTipTemplate;
        dtt.DataTipRows(1) = dataTipTextRow('r = ', sign_cc_points);
        dtt.DataTipRows(2) = dataTipTextRow('postdrug epoch', epochs);
    end
else
    h1 = plot(post_x, unsign_cc_points, 'x', 'MarkerSize', 10, 'LineWidth', 1, 'Color', [0.75 0.75 0.75]);
    hold on;
    h2 = plot(post_x, sign_cc_points, 'rx', 'MarkerSize', 10, 'LineWidth', 3);
    hold off;
    ylim([-1, 1])

    %Insert data tips
    enableDefaultInteractivity(gca);
    if ~isempty(h1)
        dtt = h1.DataTipTemplate;
        dtt.DataTipRows(1) = dataTipTextRow('r = ', unsign_cc_points);
        dtt.DataTipRows(2) = dataTipTextRow('postdrug epoch', epochs);
    end
    if ~isempty(h2)
        dtt = h2.DataTipTemplate;
        dtt.DataTipRows(1) = dataTipTextRow('r = ', sign_cc_points);
        dtt.DataTipRows(2) = dataTipTextRow('postdrug epoch', epochs);
    end
end
set(gcf, 'Position', get(0, 'Screensize'));
