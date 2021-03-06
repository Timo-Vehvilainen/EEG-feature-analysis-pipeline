%Timo Vehvilšinen November 2019

function [] = ...
    plot_artifact_percentages(predrug_artifact_epochs,postdrug_artifact_epochs, APTs, save_path)
%PLOT_ARTIFACT_PERCENTAGES 
%   Gives a visualization of which features were
%   able to be calculated for which epochs

APT_names= fieldnames(APTs);
APT_values = cellfun(@(x)(APTs.(x)),APT_names);
unique_thresholds = unique(APT_values);

%creating a gradient between green and red, with a step for each unique
%threshold
colormap_max = min([max(unique_thresholds)+0.1, 1]);
length = numel(unique_thresholds);
red = [1, 0, 0];
green = [0, 1, 0];
black = [0, 0, 0];
colors =   [linspace(green(1),red(1),length)', ...
            linspace(green(2),red(2),length)', ...
            linspace(green(3),red(3),length)'];

%Creating a colormap and tick labels for a colorbar
my_colormap = repmat(colors(1, :), floor(unique_thresholds(1)*100), 1);
threshold_labels = cell(1, length);
threshold_labels{1} = '0%';
for i=2:length
    my_colormap = [my_colormap;
                   repmat(colors(i, :), floor((unique_thresholds(i)-unique_thresholds(i-1))*100), 1)];
    features_in_range = {APT_names{APT_values == unique_thresholds(i)}};
    Str = sprintf('%s,', features_in_range{:});
    Str(end) = [];
    threshold_labels{i} = sprintf("%d%% %s", (unique_thresholds(i)*100),...
        Str);
end
my_colormap = [my_colormap;
          repmat(black, floor((colormap_max - max(unique_thresholds))*100), 1)];

predrug_artifact_prcnt = get_artifact_prcnt(predrug_artifact_epochs);
postdrug_artifact_prcnt = get_artifact_prcnt(postdrug_artifact_epochs);
           
all_subs = numel(predrug_artifact_prcnt);
f = figure;
ch_names = {'F3','F4','P3','P4','Left','Right','Frontal','Parietal'};

%Going through the subjects, doing 2 plots for each (predrug & postdrug),
%each containing a color-coded row for each channel
for sub = 1:all_subs
    
    %Predrug
     subplot(all_subs, 2, ((sub-1)*2) + 1);
     ch_no = size(predrug_artifact_prcnt{sub}, 2);
     epoch_no = size(predrug_artifact_prcnt{sub}, 3);
     predrug_image = zeros(ch_no, epoch_no);
     for epoch = 1:epoch_no
        temp = diag(squeeze(predrug_artifact_prcnt{sub}(:, :, epoch)));
        predrug_image(:, epoch) = temp(1:8);
     end
    imagesc(predrug_image, [0 colormap_max])
    ax = gca;
%     Set where ticks will be
    ax.YTick = 1:ch_no;
    ax.XTick = 1:epoch_no;
%     Set TickLabels;
    ax.YTickLabel = ch_names;
    ax.XTickLabel = split(num2str(-1*(epoch_no):-1));
    colormap(my_colormap);
    if sub == 1
        title('Artifact %-age per epoch (predrug)');
    elseif sub == all_subs
        xlabel("epoch");
    end
    
    %Postdrug
    subplot(all_subs, 2, ((sub-1)*2) + 2);
    ch_no = size(postdrug_artifact_prcnt{sub}, 2);
    epoch_no = size(postdrug_artifact_prcnt{sub}, 3);
    postdrug_image = zeros(ch_no, epoch_no);
    for epoch = 1:epoch_no
        temp = diag(squeeze(postdrug_artifact_prcnt{sub}(:, :, epoch)));
        postdrug_image(:, epoch) = temp(1:8);
    end
    imagesc(postdrug_image, [0 colormap_max])
    ax = gca;
    ax.YTick = 1:ch_no;
    ax.XTick = 1:epoch_no;
%     Set TickLabels;
    ax.YTickLabel = ch_names;
    colormap(my_colormap);
    if sub == 1
        title('Artifact %-age per epoch (postdrug)');
    elseif sub == all_subs
        xlabel("epoch");
    end
end
set(gcf, 'Position', get(0, 'Screensize'));
bottom_plot_pos = get(subplot(all_subs,2, all_subs*2),'Position');
second_bottom_plot_pos = get(subplot(all_subs,2, (all_subs-1)*2),'Position');
space_between_plots = second_bottom_plot_pos(2) - (bottom_plot_pos(2) + bottom_plot_pos(4));

%Set the colorbar. You may choose to adjust the positioning when the number
%of subjects changes
colorbar('Position', [bottom_plot_pos(1)+bottom_plot_pos(3)+0.01, ...
                        bottom_plot_pos(2), ...
                        bottom_plot_pos(4)/5, ...
                        (bottom_plot_pos(4)*all_subs + space_between_plots*(all_subs-1))], ...
        'Ticks',unique_thresholds,...
        'TickLabels',threshold_labels);

savefig(f, fullfile(save_path, 'APT.fig'));
%export_fig(f, fullfile(save_path, 'APT.pdf'), '-transparent');
%close all;
end

