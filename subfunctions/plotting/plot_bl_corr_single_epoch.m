%Timo Vehviläinen August 2019

function [] = plot_bl_corr_single_epoch(pre, post, baseline_fcn, ...
    postdrug_epoch)
%PLOT_BL_CORR_SINGLE_EPOCH
%   This function takes in epochs of pre- & postdrug EEG data for a single
% feature, calculates a baseline from the predrug epochs, calculates the
% correlation coeffiecient for a single postdrug epoch against said baseline
% (and its p-value), and finally the baseline-adjusted postdrug values
% against the baseline, and also plots the linear regression line.
%PARAMETERS:
%       pre:
%           A cell array, where each cell contains a vector or numeric
%           data, representing predrug feature values of epoched EEG data.
%       post:
%           A cell array, where each cell contains a vector or numeric
%           data, representing postdrug feature values of epoched EEG data.
%       baseline_fcn:
%           A function variable, which takes a single numeric vector as
%           it's argument. This funciton is used for calculating the
%           baseline from a single predrug feature vector.
%       postdrug_epoch:
%           a positive number, indicating which postdrug epoch is supposed
%           to be used for calculating and plotting the correlation.
%       p_threshold:
%           a value between 0 and 1, representing the significance
%           threshold chosen for the correlations. Correlations where the
%           p-value is smaller than p_threshold will be highlighted in the
%           graph.
 
all_baselines = cellfun(@(x)baseline_fcn(x), pre);
postdrug_values = cellfun(@(x)x(postdrug_epoch), post) - all_baselines;

poly = polyfit(all_baselines(~isnan(postdrug_values)),...
    postdrug_values(~isnan(postdrug_values)),1);

% Evaluate the fitted polynomial and plot:
outliers = isoutlier(postdrug_values);
f = polyval(poly, all_baselines);
plot(all_baselines,f,'-', 'color',[0.9100    0.4100    0.1700]);
hold on;
plot(all_baselines(~outliers), postdrug_values(~outliers), 'bo');
plot(all_baselines(outliers), postdrug_values(outliers), 'rx');

%Type of correlation to calculate
correlation_type = 'Spearman';
corr_label = {"Spearman's \rho"};

[CC, P] = corr([postdrug_values, all_baselines], 'rows','pairwise', 'Type', correlation_type);
p = P(1, 2);
legend(sprintf('%s = %.3f\np = %.3f', corr_label{1},  CC(1, 2), p));
