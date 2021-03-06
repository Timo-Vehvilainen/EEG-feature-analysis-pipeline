%Timo Vehvilšinen August 2019

function [] = plot_bl_corr(pre, post, epoch_length_seconds, epoch_overlap_seconds, ...
    max_epochs, baseline_fcn, p_threshold)
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
%       max_epochs:
%           A positive integer for limiting the number of epochs that are plotted
%       baseline_fcn:
%           A function variable, which takes a single numeric vector as
%           it's argument. This funciton is used for calculating the
%           baseline from a single predrug feature vector.
%       p_threshold:
%           a value between 0 and 1, representing the significance
%           threshold chosen for the correlations. Correlations where the
%           p-value is smaller than p_threshold will be highlighted in the
%           graph.

max_postdrug_epochs = min(max_epochs, max(cellfun(@numel, post)));
post_x = (epoch_length_seconds/2): ...
    (epoch_length_seconds - epoch_overlap_seconds): ...
    ((epoch_length_seconds*(max_epochs)/2));

all_baselines = cellfun(@(x)baseline_fcn(x), pre);
[ccs, ps] = deal(zeros(1, max_postdrug_epochs));
for epoch = 1:max_postdrug_epochs

    post_values = cellfun(@(x)x(epoch), post) - all_baselines;
    %remove the single highest and lowest values to lessen the impact of
    %outliers
    post_values(post_values==max(post_values)) = NaN;
    post_values(post_values==min(post_values)) = NaN;
    [CC, P] = corr([post_values, all_baselines], ...
        'rows','pairwise');
    ccs(epoch) = CC(1, 2);
    ps(epoch) = P(1, 2);
end

p_significant = (ps < p_threshold);

%Plot a horizontal line at the chosen significance level
epochs = 1:max_postdrug_epochs;
yyaxis right
plot([0, max(post_x)], [p_threshold, p_threshold], '--');
hold on;

%Plot the p-values (differently depending on if they are below the
%significance level)
sign_points = ps;
sign_points(~p_significant) = NaN;
unsign_points = ps;
unsign_points(p_significant) = NaN;

h1 = plot(post_x, unsign_points, 'rx', 'LineWidth', 1, 'MarkerSize', 10);
hold on
h2 = plot(post_x, sign_points, 'ro', 'LineWidth', 3, 'MarkerSize', 10);
hold off;

%insert data tips
enableDefaultInteractivity(gca);
if ~isempty(h1)
    dtt = h1.DataTipTemplate;
    dtt.DataTipRows(1) = dataTipTextRow('p =', unsign_points);
    dtt.DataTipRows(2) = dataTipTextRow('postdrug epoch', epochs);
end
if ~isempty(h2)
    dtt = h2.DataTipTemplate;
    dtt.DataTipRows(1) = dataTipTextRow('p =', sign_points);
    dtt.DataTipRows(2) = dataTipTextRow('postdrug epoch', epochs);
end

%Plot the correlation coefficients (differently depending on if they have a
%significant p-value or not).
yyaxis left
sign_points = ccs;
sign_points(~p_significant) = NaN;
unsign_points = ccs;
unsign_points(p_significant) = NaN;

h3 = plot(post_x, unsign_points, 'bx', 'MarkerSize', 10, 'LineWidth', 1);
hold on;
h4 = plot(post_x, sign_points, 'bo', 'MarkerSize', 10, 'LineWidth', 3);
hold off;
set(gca, 'YDir','reverse')

%Insert data tips
enableDefaultInteractivity(gca);
if ~isempty(h3)
    dtt = h3.DataTipTemplate;
    dtt.DataTipRows(1) = dataTipTextRow('cc = ', unsign_points);
    dtt.DataTipRows(2) = dataTipTextRow('postdrug epoch', epochs);
end
if ~isempty(h4)
    dtt = h4.DataTipTemplate;
    dtt.DataTipRows(1) = dataTipTextRow('cc = ', sign_points);
    dtt.DataTipRows(2) = dataTipTextRow('postdrug epoch', epochs);
end
set(gcf, 'Position', get(0, 'Screensize'));
