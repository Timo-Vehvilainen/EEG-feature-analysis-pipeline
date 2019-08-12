%Timo Vehviläinen August 2019

function [features] = get_features(epoched_data, artifact_epochs, varargin)
%GET_FEATURES Calculates computational features from epoched EEG data
%   This function takes in epoched multi-channel EEG data from multiple 
%   subjects and calculates a set of predefined features for that data. The
%   features calculated are the following (and the amount of features
%   calculated for each subject):
%                 Activation synchrony Index (ASI),   
%                       3D array: [channels * channels * epochs]
%                 Weighted Phase Lag Index (wPLI),
%                       4D array: [channels * channels * wPLI_frequency_bands * epochs]
%                 Nestedness coefficient (NC),      
%                       3D array: [channels * NC_frequency_bands * epochs]
%                 amplitude-integrated EEG (aEEG),
%                       - subfeatures:  
%                           mean, 
%                           interquartile range (iqr)
%                               2D array: [channels * epochs] for each
%                               subfeature
%                 range-EEG (rEEG),
%                       - subfeatures:
%                           mean,
%                           interquartile range (iqr),
%                           lower 5th percentile (li)
%                               2D array: [channels * epochs] for each
%                               subfeature
%                 Power Spectral Density (PSD),  
%                       3D array: [channels * PSD_frequency_bands * epochs]
%                 Cross Power Spectral Density (cPSD),     
%                       4D array: [channels * channels * PSD_frequency_bands * epochs]
%                 Multifractal Detrended Fluctuation analysis (MFDFA),
%                       - subfeatures:
%                           width,
%                           height,
%                           tail,
%                           peak
%                               2D array: [channels * epochs] for each subfeature
%                 Suppression curve (SC)
%                       2D array: [channels * epochs]
% PARAMETERS:
%           epoched_data:   a cell vector of epoched EEG data, provided by
%               the epoch_data()-function. Each cell element is a 3D-array of
%               the form [time_samples * channels * epochs] and represents the
%               epoched data for a single EEG subject/patient.
%           artifact_epochs:    a cell vector similar to epoched_data, but
%               one where the data is binary artifact mask denoting the periods
%               of the data that have been deemed to contain artifact or are
%               missing data. 1s represent missing data, 0s represent good
%               data. For each epoch, the percentage of compromised samples
%               is calculated and compared against an internally defined
%               APT value (view variable artifact_percentage_thresholds below).
%               If the detected percentage for a given epoch 
%               is greater than the predefined value for that feature, the
%               calculation of that feature is omitted for that epoch.
%
%OPTIONAL PARAMETERS:
%           epoch_length_seconds: The length of the epochs in seconds.
%           Defaults to 120 seconds.
%
%           sampling_rate: The sampling rate frequency of the data
%           provided. Defaults to 250 Hz.
%
%OUTPUT is a struct where each field is a cell vector named after a
%feature. Each cell element represents a subject and is of the form specified
%above for the given feature.

%parse input
p = inputParser;
addParameter(p, 'epoch_length_seconds', 120, @isnumeric);
addParameter(p, 'sampling_rate', 250, @isnumeric);
addParameter(p, 'MFDFA_scale_min', 25, @isnumeric);
addParameter(p, 'MFDFA_scale_max', 1875, @isnumeric);
addParameter(p, 'MFDFA_scale_res', 19, @isnumeric);
addParameter(p, 'MFDFA_q', -5:0.5:5, @isnumeric);
addParameter(p, 'MFDFA_m', 1, @isnumeric);

parse(p, varargin{:});
S = struct2cell(p.Results); 
[epoch_length_seconds, MFDFA_m, MFDFA_q, MFDFA_scale_max, MFDFA_scale_min, MFDFA_scale_res, ...
    sampling_rate] = deal(S{:});

all_subs = numel(epoched_data);
ch_no = 8;

%Get the percentage of artifacts per epoch
artifact_prcnts = get_artifact_prcnt(artifact_epochs);

%Define artifact percentage thresholds (APT)
artifact_percentage_thresholds = ...
                    struct( 'ASI',     0.2, ...
                            'wPLI',     0.1, ...
                            'NC',       0.1, ...
                            'aEEG',     0.5, ...
                            'rEEG',     0.5, ...
                            'PSD',      0.5, ...
                            'cPSD',     0.5, ...
                            'MFDFA',    0.01, ...
                            'SC',       0.5);

%define frequency bands for wPLI, NC and PSD
NC_highcomps = [3, 8;
                8, 15;
                15, 30];
NC_band_no = size(NC_highcomps, 1);
wPLI_bands = [0.4, 3;
            4, 8;
            8, 13;
            13, 22];
wPLI_band_no = size(wPLI_bands, 1);
PSD_bands =    [1, 3;
                3, 8;
                8, 15;
                15, 30];
PSD_band_no = size(PSD_bands, 1);

%define settings for calculating line length
window_size_seconds = 1;
overlap_seconds = 0.12;
window_size_samples = floor(window_size_seconds * sampling_rate);
overlap_samples = floor(overlap_seconds * sampling_rate);
window_interval = window_size_samples - overlap_samples;
epoch_length_samples = epoch_length_seconds * sampling_rate;
Ll_size = floor((epoch_length_samples - window_interval) / overlap_samples);

%initiate features
[ASI, wPLI, NC, aEEG, rEEG, PSD, cPSD, MFDFA, SC] = deal(cell(all_subs, 1));
MFDFA_exponents = linspace(log2(MFDFA_scale_min),log2(MFDFA_scale_max),MFDFA_scale_res);
MFDFA_scales = round(2.^MFDFA_exponents);

%loop through subjects
for sub = 1:all_subs
    multiWaitbar(sprintf('Subjects (%d total)', all_subs), 'Value', sub/all_subs, 'Color', 'g');
    if isnan(epoched_data{sub})
        [ASI{sub}, wPLI{sub}, NC{sub}, aEEG{sub}, rEEG{sub}, ...
            PSD{sub}, cPSD{sub}, MFDFA{sub}, SC(sub)] = deal(NaN);
        continue;
    end
    
    %Initiate features for one subject
    epoch_no = size(epoched_data{sub}, 3);
    ASI{sub} = zeros(ch_no, ch_no, epoch_no);
    cPSD{sub} = zeros(ch_no, ch_no, PSD_band_no, epoch_no);
    wPLI{sub} = zeros(ch_no, ch_no, wPLI_band_no, epoch_no);
    [aEEG{sub}.mean, aEEG{sub}.iqr, ...
        rEEG{sub}.mean, rEEG{sub}.iqr, rEEG{sub}.li, ...
        MFDFA{sub}.width, MFDFA{sub}.height, ...
        MFDFA{sub}.peak, MFDFA{sub}.tail] = deal(zeros(8, epoch_no));
    PSD{sub} = zeros(ch_no, PSD_band_no, epoch_no);
    NC{sub} = zeros(ch_no, NC_band_no, epoch_no);
    SC{sub} = zeros(1, epoch_no);
    
    %loop through all the epochs for one subject
    line_lengths = zeros(Ll_size, 4);
    for epoch = 1:epoch_no
        multiWaitbar('Epochs for this subject...', 'Value', epoch/epoch_no, 'Color', 'g');
        active_epoch = squeeze(epoched_data{sub}(:, :, epoch));
        artifact_epoch = squeeze(artifact_epochs{sub}(:, :, epoch));
        
        %loop through the channels for one subject
        for ch1 = 1:ch_no
            for ch2 = 1:ch_no
                %calculate values for cross-channel features (ASI, wPLI & cPSD)
                artifact_prcnt = artifact_prcnts{sub}(ch1, ch2, epoch);
                
                %the artifact masks from both channels are combined before
                %application
                combined_artifact_mask = or(artifact_epoch(:, ch1), artifact_epoch(:, ch2));
                masked_epoch_ch1 = active_epoch(:, ch1);
                masked_epoch_ch1 = masked_epoch_ch1(~combined_artifact_mask);
                masked_epoch_ch2 = active_epoch(:, ch2);
                masked_epoch_ch2 = masked_epoch_ch2(~combined_artifact_mask); 
                
                %wPLI
                if artifact_prcnt <= artifact_percentage_thresholds.wPLI
                    %looping through the wPLI frequency bands
                    for wpli_hz = 1:wPLI_band_no
                        bandpassed_epoch_ch1 = my_bandpass(masked_epoch_ch1, wPLI_bands(wpli_hz, :), sampling_rate);
                        bandpassed_epoch_ch2 = my_bandpass(masked_epoch_ch2, wPLI_bands(wpli_hz, :), sampling_rate);
                        wPLI{sub}(ch1, ch2, wpli_hz, epoch) = get_wPLI(hilbert(bandpassed_epoch_ch1), ...
                                                        hilbert(bandpassed_epoch_ch2), 1, wPLI_bands(wpli_hz, :), sampling_rate);
                    end
                else
                        wPLI{sub}(ch1, ch2, :, epoch) = NaN;
                end
                
                %ASI
                if artifact_prcnt <= artifact_percentage_thresholds.ASI
                    %bandpass the epoch between [1.5 20] Hz for ASI
                    asi_bandpassed_epoch_ch1 = my_bandpass(masked_epoch_ch1, [1.5 20], sampling_rate);
                    asi_bandpassed_epoch_ch2 = my_bandpass(masked_epoch_ch2, [1.5 20], sampling_rate);
                    ASI{sub}(ch1, ch2, epoch) = getASI([asi_bandpassed_epoch_ch1, ...
                                                    asi_bandpassed_epoch_ch2], sampling_rate);
                else
                    ASI{sub}(ch1, ch2, epoch) = NaN;
                end
                
                %cPSD
                if artifact_prcnt < artifact_percentage_thresholds.cPSD
                    %loop through the PSD frequency bands
                    for cPSD_band = 1:PSD_band_no
                        [b,a] = butter(5, PSD_bands(cPSD_band, 2)/sampling_rate, 'low'); 
                        filtered_epoch_ch11 = filtfilt(b,a,masked_epoch_ch1); 
                        [b,a] = butter(5, PSD_bands(cPSD_band, 1)/sampling_rate, 'high');
                        filtered_epoch_ch11 = filtfilt(b, a, filtered_epoch_ch11);

                        [b,a] = butter(5, PSD_bands(cPSD_band, 2)/sampling_rate, 'low'); 
                        filtered_epoch_ch2 = filtfilt(b,a,masked_epoch_ch2); 
                        [b,a] = butter(5, PSD_bands(cPSD_band, 1)/sampling_rate, 'high');
                        filtered_epoch_ch2 = filtfilt(b, a, filtered_epoch_ch2);

                        cPSD{sub}(ch1, ch2, cPSD_band, epoch) = abs(mean(cpsd(filtered_epoch_ch11, filtered_epoch_ch2, ...
                                                10*sampling_rate,0,20*sampling_rate)));
                
                    end
                else
                        cPSD{sub}(ch1, ch2, :, epoch) = NaN;
                end
            end
            artifact_prcnt = artifact_prcnts{sub}(ch1, ch1, epoch);
            masked_epoch_ch1 = active_epoch(:, ch1);
            masked_epoch_ch1 = masked_epoch_ch1(~artifact_epoch(:, ch1));
            
            %calculating non cross-channel features
            
            if artifact_prcnt <= artifact_percentage_thresholds.MFDFA
                %MFDFA
                [~, ~, hq, Dq, ~] = MFDFA1(masked_epoch_ch1, MFDFA_scales, MFDFA_q, MFDFA_m, 0);
                [Dq_max, Dq_max_i] = max(Dq);
                [hq_max, hq_max_i] = max(hq);
                [~, hq_min_i] = min(hq);
                MFDFA{sub}.width(ch1, epoch) = hq_max - min(hq);
                MFDFA{sub}.height(ch1, epoch) = Dq_max - min(Dq);
                MFDFA{sub}.peak(ch1, epoch) = hq(Dq_max_i);
                MFDFA{sub}.tail(ch1, epoch) = Dq(hq_min_i) - Dq(hq_max_i);
            else
                MFDFA{sub}.width(ch1, epoch) = NaN;
                MFDFA{sub}.height(ch1, epoch) = NaN;
                MFDFA{sub}.peak(ch1, epoch) = NaN;
                MFDFA{sub}.tail(ch1, epoch) = NaN;
            end
            
            % Nestedness
            if artifact_prcnt <= artifact_percentage_thresholds.NC
                %loop through NC frequency bands
                for NC_band = 1:NC_band_no
                    NC{sub}(ch1, NC_band, epoch) = get_nestedness(masked_epoch_ch1, NC_highcomps(NC_band, :), sampling_rate);
                end
            else
                    NC{sub}(ch1, :, epoch) = NaN;
            end
            
            %aEEG
            if artifact_prcnt <= artifact_percentage_thresholds.aEEG
                [~, aeeg] = my_aeeg(masked_epoch_ch1, sampling_rate, 0);
                aEEG{sub}.mean(ch1, epoch) = mean(aeeg);
                aEEG{sub}.iqr(ch1, epoch) = iqr(aeeg);
            else
                aEEG{sub}.mean(ch1, epoch) = NaN;
                aEEG{sub}.iqr(ch1, epoch) = NaN;
            end
            
            %rEEG
            if artifact_prcnt <= artifact_percentage_thresholds.rEEG
                [~, reeg] = estimate_rEEG(masked_epoch_ch1, sampling_rate);
                rEEG{sub}.mean(ch1, epoch) = mean(reeg);
                rEEG{sub}.iqr(ch1, epoch) = iqr(reeg);
                rEEG{sub}.li(ch1, epoch) = prctile(reeg, 5);
            else
                rEEG{sub}.mean(ch1, epoch) = NaN;
                rEEG{sub}.iqr(ch1, epoch) = NaN;
                rEEG{sub}.li(ch1, epoch) = NaN;
            end
            
            %PSD
            if artifact_prcnt <= artifact_percentage_thresholds.PSD
                %loop through PSD frequency bands
                for PSD_band = 1:PSD_band_no
                    [b,a] = butter(5, PSD_bands(PSD_band, 2)/sampling_rate, 'low'); 
                    filtered_epoch = filtfilt(b,a,masked_epoch_ch1); 
                    [b,a] = butter(5, PSD_bands(PSD_band, 1)/sampling_rate, 'high');
                    filtered_epoch = filtfilt(b, a, filtered_epoch);
                    PSD{sub}(ch1, PSD_band, epoch) = mean(pwelch(filtered_epoch,...
                                        10*sampling_rate,[],sampling_rate*20));
                end
            else
                PSD{sub}(ch1, :, epoch) = NaN;
            end
            
            % Line length (for calculating suppression curve)
            if artifact_prcnt <= artifact_percentage_thresholds.SC 
                %line lenght is only calculated for monopolar channels
                if ch1 < 5
                    line_length_epoch = [masked_epoch_ch1;
                                        NaN(epoch_length_samples - numel(masked_epoch_ch1), 1)];
                    line_lengths(:, ch1) = get_line_length(line_length_epoch, 'window_size_samples', window_size_samples, ...
                                                        'overlap_samples', overlap_samples); 
                end
            else
                line_lengths(:, ch1) = NaN;
            end
        end
        
        %Extract the median of the normalized line lengths for the monopolar channels, 
        %and calculate the suppression curve 
        median_curve = nanmedian(line_lengths(:, 1:4), 2)';
        median_curve(median_curve == 0) = [];
        suppr_curve_value = 1 - (nanmedian(median_curve) / nanmean(median_curve));
        SC{sub}(epoch) = suppr_curve_value;
    end
    multiWaitbar('Epochs for this subject...', 'Close');
end
multiWaitbar('CLOSEALL');

%compile all the collected features into a struct
features = struct('ASI', ASI, 'wPLI', wPLI, 'NC', NC, 'aEEG', ...
                    aEEG, 'rEEG', rEEG, 'PSD', PSD, 'cPSD', cPSD, ...
                    'MFDFA', MFDFA, 'SC', SC);
% close(f)
end

