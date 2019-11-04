%Timo Vehviläinen August 2019

function [data_table] = plot_bl_timeline(pre, post, epoch_length_seconds, epoch_overlap_seconds, ...
    show_epochs, baseline_fcn, baseline_mode)
%PLOT_TIMELINE 
%   This function takes in epochs of pre- & postdrug EEG data for a single
%   feature, calculates a baseline from the predrug epochs, and depending
%   on the given baseline mode, plots the ratio or difference of the EEG
%   time trends against the baseline. It plots them in a graph, and 
%   additionally plots a trimmed mean signal as well (the highest and the lowest
%   value are ignored when calculating the mean trend).
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
%   Returns a table of variables. The first column variable holds the
%   x-axis values for the plot, and the rest hold the y-data for the
%   plotted lines, including the mean.

max_predrug_epochs = min(show_epochs(1), max(cellfun(@numel, pre)));
max_postdrug_epochs = min(show_epochs(2), max(cellfun(@numel, post)));
pre_values_all = NaN(numel(pre), max_predrug_epochs);
post_values_all = NaN(numel(post), max_postdrug_epochs);
x_limits = [-1*epoch_length_seconds * (max_predrug_epochs)/2, epoch_length_seconds * (max_postdrug_epochs-1)/2];

pre_x = -1 *((epoch_length_seconds/2):...
            epoch_overlap_seconds:...
            ((epoch_length_seconds*max_predrug_epochs)/2));

post_x = (epoch_length_seconds/2):...
            epoch_overlap_seconds:...
            ((epoch_length_seconds*max_postdrug_epochs)/2);        
        
x_axis = [fliplr(pre_x)'; post_x'];
data_table = table(x_axis, 'VariableNames', "X_axis");

for sub = 1:numel(pre)
    subject_name = sprintf("Subject_%02d", sub);
    if ~isempty(pre{sub})
        pre_values = squeeze(pre{sub})';
        pre_values_no = numel(pre_values);
        pre_values = pre_values(max(pre_values_no-max_predrug_epochs+1, 1):pre_values_no);
        baseline = baseline_fcn(pre_values);
        if strcmp(baseline_mode, 'relative')
            pre_values = pre_values ./ baseline;
            %change unit to percentages:
            pre_values = pre_values * 100;
        else
            pre_values = pre_values - baseline;
        end
    else
        continue;
    end
    pre_values_all(sub, :) = [NaN(1, (size(pre_values_all, 2) - numel(pre_values))), pre_values];
    if ~isempty(post{sub})
        if strcmp(baseline_mode, 'relative')
            post_values = squeeze(post{sub}(1:min(max_postdrug_epochs, end)))'./baseline;
            %change value to percentages:
            post_values = post_values * 100;
        else
            post_values = squeeze(post{sub}(1:min(max_postdrug_epochs, end)))' - baseline;
        end
    else
        post_values = [];
    end
    post_values_all(sub, :) = [post_values, NaN(1, (size(post_values_all, 2) - size(post_values, 2)))];
    plot_value_no = numel([pre_x, post_x]);
    l = plot([fliplr(pre_x), post_x], [pre_values, post_values], 'color', [0.6, 0.6, 0.6], 'LineWidth', 1);
    hold on;
    xlim(x_limits);
    
    data_table = addvars(data_table, [pre_values'; post_values'], 'NewVariableNames', subject_name);
    
    if strcmp(baseline_mode, 'relative')
        ytickformat('percentage')
    end
    enableDefaultInteractivity(gca);
    if ~isempty(l)
        dtt = l.DataTipTemplate;
        dtt.DataTipRows(1) = dataTipTextRow('epoch', [(-1*(numel(pre_x)):numel(post_x))+1]);
        dtt.DataTipRows(2) = dataTipTextRow('feature value', [pre_values, post_values]);
        dtt.DataTipRows(end+1) = dataTipTextRow('Subject', ones(plot_value_no, 1)*sub);
    end
end

%Plot trimmed mean trend
trim_percentage = 200/numel(pre);
pre_mean = trimmean(pre_values_all, trim_percentage);
post_mean = trimmean(post_values_all, trim_percentage);
pre_mean_x = -1 *((epoch_length_seconds/2):epoch_overlap_seconds:(epoch_length_seconds*(numel(pre_mean)))/2);
post_mean_x = (epoch_length_seconds/2):epoch_overlap_seconds:((epoch_length_seconds*(numel(post_mean)))/2);

data_table = addvars(data_table, [pre_mean'; post_mean'], 'NewVariableNames', "TrimmedMean");

plot([fliplr(pre_mean_x), post_mean_x], [pre_mean, post_mean], 'LineWidth', 3, 'color', 'blue');
set(gcf, 'Position', get(0, 'Screensize'));

vline(0, [0, 0, 0]);
end

