function [baseline] = baseline_median5(predrug_feature_epochs)
%BASELINE_MEDIAN5 Summary of this function goes here
%   Detailed explanation goes here
    baseline = nanmedian(predrug_feature_epochs(max(1, (end-5)):end));
end

