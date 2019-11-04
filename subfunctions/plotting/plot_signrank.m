function [data_table] = plot_signrank(pre, post, ...
    epoch_length_seconds, epoch_overlap_seconds, ...
    show_epochs, baseline_fcn, p_threshold, hide_p_values)
%PLOT_SIGNRANK
%   This function takes in epochs of pre- & postdrug EEG data for a single
%   feature, calculates the Wilcoxon signed-rank test for all the postdrug
%   epochs against a chosen baseline. It plots them in a timeline graph with
%   their p-values, which can optionally be signified with just highlights 
%   using the hide_p_values variable
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
%           are to be plotted from the data.
%       baseline_fcn:
%           A function variable, which takes a single numeric vector as
%           it's argument. This funciton is used for calculating the
%           baseline from a single predrug feature vector.
%       baseline_mode:
%           a string or character vector, which is either 'relative' or
%           'absolute'. If the mode is set to relative, then the ratio of
%           (time-trend/baseline) is plotted. If set to absolute, then the
%           difference of (time-trend - baseline) is plotted.
%OUTPUT
%   Return a table of column variables, where the first variable is the
%   x-axis values used to plot the signranks. The second variable holds the
%   calculated median values for each post-drug epoch, and the third
%   variable holds their p-values


max_postdrug_epochs = min(show_epochs(2), max(cellfun(@numel, post)));
post_x = (epoch_length_seconds/2): ...
    (epoch_length_seconds - epoch_overlap_seconds): ...
    ((epoch_length_seconds*(max_postdrug_epochs)/2));

data_table = table(post_x', 'VariableNames', "X_axis");
all_baselines = cellfun(@(x)baseline_fcn(x), pre);
[medians, ps] = deal(zeros(1, max_postdrug_epochs));

for epoch = 1:max_postdrug_epochs
    post_values = cellfun(@(x)x(epoch), post);
    medians(epoch) = nanmedian(post_values - all_baselines);
    [ps(epoch), H] = signrank(all_baselines, post_values, 'alpha', p_threshold);
end

data_table = addvars(data_table, medians', ps', 'NewVariableNames', {'Medians', 'P_values'});

p_significant = (ps < p_threshold);

%Plot a horizontal line at the chosen significance level
epochs = 1:max_postdrug_epochs;

%separate the medians and p-vals into significant and unsignificant sets
sign_p_points = ps;
sign_p_points(~p_significant) = NaN;
unsign_p_points = ps;
unsign_p_points(p_significant) = NaN;

sign_median_points = medians;
sign_median_points(~p_significant) = NaN;
unsign_median_points = medians;
unsign_median_points(p_significant) = NaN;

%plot differently depending on the hide_p_values variable
%if false, p-values will be plotted on the same graph with a second y-axis
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

    plot([0, max(post_x)], [0, 0], 'b--');
    hold on;
    h3 = plot(post_x, unsign_median_points, 'bx', 'MarkerSize', 10, 'LineWidth', 1);
    h4 = plot(post_x, sign_median_points, 'bo', 'MarkerSize', 10, 'LineWidth', 3);
    hold off;

    %Insert data tips
    enableDefaultInteractivity(gca);
    if ~isempty(h3)
        dtt = h3.DataTipTemplate;
        dtt.DataTipRows(1) = dataTipTextRow('median = ', unsign_median_points);
        dtt.DataTipRows(2) = dataTipTextRow('postdrug epoch', epochs);
    end
    if ~isempty(h4)
        dtt = h4.DataTipTemplate;
        dtt.DataTipRows(1) = dataTipTextRow('median = ', sign_median_points);
        dtt.DataTipRows(2) = dataTipTextRow('postdrug epoch', epochs);
    end
%if true, plot only the medians and highlight based on significance
else
    plot([0, max(post_x)], [0, 0], '--', 'Color', [0.75 0.75 0.75]);
    hold on;
    h1 = plot(post_x, unsign_median_points, 'x', 'MarkerSize', 10, 'LineWidth', 1, 'Color', [0.75 0.75 0.75]);
    h2 = plot(post_x, sign_median_points, 'rx', 'MarkerSize', 10, 'LineWidth', 3);
    hold off;

    %Insert data tips
    enableDefaultInteractivity(gca);
    if ~isempty(h1)
        dtt = h1.DataTipTemplate;
        dtt.DataTipRows(1) = dataTipTextRow('median = ', unsign_median_points);
        dtt.DataTipRows(2) = dataTipTextRow('postdrug epoch', epochs);
    end
    if ~isempty(h2)
        dtt = h2.DataTipTemplate;
        dtt.DataTipRows(1) = dataTipTextRow('median = ', sign_median_points);
        dtt.DataTipRows(2) = dataTipTextRow('postdrug epoch', epochs);
    end
end
set(gcf, 'Position', get(0, 'Screensize'));
end

