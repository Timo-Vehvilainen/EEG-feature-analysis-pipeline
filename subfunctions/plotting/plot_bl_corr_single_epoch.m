%Timo Vehviläinen August 2019

function [] = plot_bl_corr_single_epoch(pre, post, baseline_fcn, ...
    postdrug_epoch, p_threshold)
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
postdrug_values = cellfun(@(x)x(postdrug_epoch), post);

poly = polyfit(all_baselines,postdrug_values,1);
% Evaluate the fitted polynomial and plot:

f = polyval(poly, all_baselines);
figure;
plot(all_baselines,f,'-');
[CC, P] = corr([postdrug_values, all_baselines], 'rows','complete');
p = P(1, 2);
if (p > p_threshold) || (isnan(p)) 
    legend(sprintf('cc = %.3f', CC(1, 2)), sprintf('p = %.3f', P(1, 2)));
else
    legend(sprintf('cc = %.3f', CC(1, 2)), ['{\color{red}p = ' num2str(P(1, 2), '%.3f}')]);
end