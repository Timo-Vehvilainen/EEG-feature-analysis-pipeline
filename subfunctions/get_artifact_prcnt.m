%Timo Vehviläinen August 2019

function [artifact_prcnt] = get_artifact_prcnt(epoched_artifact_mask)
%GET_ARTIFACT_PRCNT 
%Calculates the percentage of samples marked as
%'artifact' or 'missing data'. for each epoch in an epoched artifact mask,
%and returns a mapped version of that mask, where each epoch is replaced by
%its corresponding percentage value in the interval [0 1].

subs = numel(epoched_artifact_mask);
artifact_prcnt = cell(subs, 1);
for sub = 1:subs
    ch_no = size(epoched_artifact_mask{sub}, 2);
    epoch_no = size(epoched_artifact_mask{sub}, 3);
    %percentage value at position (i, j) indicates the combined artifact
    %percentage of those two channels
    artifact_prcnt{sub} = zeros(ch_no, ch_no, epoch_no);
    for epoch = 1:epoch_no
        artifact_epoch = squeeze(epoched_artifact_mask{sub}(:, :, epoch));
        epoch_length_samples = size(artifact_epoch, 1);
        for ch1 = 1:ch_no
            for ch2 = 1:ch_no
                combined_channel_mask = or(artifact_epoch(:, ch1), artifact_epoch(:, ch2));
                artifact_prcnt{sub}(ch1, ch2, epoch) = numel(find(combined_channel_mask)) / epoch_length_samples;
            end
        end 
    end
end
end

