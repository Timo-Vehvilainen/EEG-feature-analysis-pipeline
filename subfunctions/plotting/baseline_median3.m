function [baseline] = baseline_median3(predrug_feature_epochs)
%BASELINE_MEDIAN3 Summary of this function goes here
%   Detailed explanation goes here
baseline = median(predrug_feature_epochs(max(1, (end-3)):end));
end

