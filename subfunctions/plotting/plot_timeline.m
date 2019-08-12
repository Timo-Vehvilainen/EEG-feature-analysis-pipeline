%Timo Vehviläinen August 2019

function [] = plot_timeline(pre, post, epoch_length_seconds, epoch_overlap_seconds, max_epochs)
%PLOT_TIMELINE 
%   This function takes in epochs of pre- & postdrug EEG data for a single
%   feature, and plots the time trend of that feature data for multiple
%   subjects. It additionally plots a mean trend.
%PARAMETERS
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
%           A positive integer for limiting the number of epochs that are
%           plotted.

max_predrug_epochs = min(max_epochs, max(cellfun(@numel, pre)));
max_postdrug_epochs = min(max_epochs, max(cellfun(@numel, post)));
pre_values_all = NaN(numel(pre), max_predrug_epochs);
post_values_all = NaN(numel(post), max_postdrug_epochs);
x_limits = [-1*epoch_length_seconds * (max_predrug_epochs)/2, ...
    epoch_length_seconds * (max_postdrug_epochs-1)/2];

for sub = 1:numel(pre)
        %extract the predrug values
        if ~isempty(pre{sub})
            pre_values = squeeze(pre{sub})';
            pre_values = pre_values(max(end-max_epochs+1, 1):end);
        else
            pre_values = [];
        end
        pre_values_all(sub, :) = ...
            [NaN(1, (size(pre_values_all, 2) - numel(pre_values))), pre_values];
        pre_x = -1 *((epoch_length_seconds/2):...
            epoch_overlap_seconds:...
            ((epoch_length_seconds*(numel(pre_values)))/2));
        
        %extract postdrug values
        if ~isempty(post{sub})
            post_values = squeeze(post{sub})';
            post_values = post_values(1:min(end, max_epochs));
        else
            post_values = [];
        end
        post_values_all(sub, :) = ...
            [post_values, NaN(1, (size(post_values_all, 2) - size(post_values, 2)))];
        post_x = (epoch_length_seconds/2):...
            epoch_overlap_seconds:...
            ((epoch_length_seconds*(numel(post_values)))/2);
        plot_value_no = numel([pre_x, post_x]);
        
        %Plot together on an x-axis centered on the drug administration
        %point
        l = plot([fliplr(pre_x), post_x], [pre_values, post_values], 'color', [0.6, 0.6, 0.6], 'LineWidth', 1);
        hold on;
        xlim(x_limits);
        
        %Add data tips to graphed data
        enableDefaultInteractivity(gca);
        if ~isempty(l)
            dtt = l.DataTipTemplate;
            dtt.DataTipRows(1) = dataTipTextRow('epoch', [(-1*(numel(pre_x)):numel(post_x))+1]);
            dtt.DataTipRows(2) = dataTipTextRow('feature value', [pre_values, post_values]);
            dtt.DataTipRows(end+1) = dataTipTextRow('Subject', ones(plot_value_no, 1)*sub);
        end
end

%calculating and plotting the mean trend
pre_mean = nanmean(pre_values_all);
post_mean = nanmean(post_values_all);
pre_mean_x = -1 *((epoch_length_seconds/2):...
    epoch_overlap_seconds:...
    (epoch_length_seconds*(numel(pre_mean)))/2);
post_mean_x = (epoch_length_seconds/2):...
    epoch_overlap_seconds:...
    ((epoch_length_seconds*(numel(post_mean)))/2);

plot([fliplr(pre_mean_x), post_mean_x], [pre_mean, post_mean], 'LineWidth', 3, 'color', 'red');
set(gcf, 'Position', get(0, 'Screensize'));

vline(0, [0, 0, 0]);
end

