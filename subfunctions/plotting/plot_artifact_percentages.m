%Timo Vehviläinen August 2019

function [] = ...
    plot_artifact_percentages(predrug_artifact_epochs,postdrug_artifact_epochs, save_path)
%PLOT_ARTIFACT_PERCENTAGES Gives a visualization of which features were
%able to be calculated for which epochs
%   

predrug_artifact_prcnt = get_artifact_prcnt(predrug_artifact_epochs);
postdrug_artifact_prcnt = get_artifact_prcnt(postdrug_artifact_epochs);

my_colormap = [repmat([0 1 0], 1, 1);
               repmat([127/255 1 212/255], 9, 1);
               repmat([0 0 1], 10, 1); 
               repmat([1 0 0], 30, 1);
               repmat([0 0 0], 10, 1)];
           
all_subs = numel(predrug_artifact_prcnt);
f = figure;

for sub = 1:all_subs
     subplot(all_subs, 2, ((sub-1)*2) + 1);
     ch_no = size(predrug_artifact_prcnt{sub}, 1);
     epoch_no = size(predrug_artifact_prcnt{sub}, 3);
     predrug_image = zeros(8, epoch_no);
     for epoch = 1:epoch_no
        temp = diag(squeeze(predrug_artifact_prcnt{sub}(:, :, epoch)));
        predrug_image(:, epoch) = temp(1:8);
     end
    imagesc(predrug_image, [0 0.6])
    ax = gca;
%     Set where ticks will be
    ax.YTick = 1:8;
    ax.XTick = 1:epoch_no;
%     Set TickLabels;
    ax.YTickLabel = {'F3','F4','P3','P4','Left','Right','Frontal','Parietal'};
    colormap(my_colormap);
    if sub == 1
        title('Artifact %-age per epoch, (Subjects 1-4, predrug)');
    elseif sub == all_subs
        xlabel("epoch");
    end
    
    
    subplot(all_subs, 2, ((sub-1)*2) + 2);
    ch_no = size(postdrug_artifact_prcnt{sub}, 1);
    epoch_no = size(postdrug_artifact_prcnt{sub}, 3);
    postdrug_image = zeros(8, epoch_no);
    for epoch = 1:epoch_no
        temp = diag(squeeze(postdrug_artifact_prcnt{sub}(:, :, epoch)));
        postdrug_image(:, epoch) = temp(1:8);
    end
    imagesc(predrug_image, [0 0.6])
    ax = gca;
    ax.YTick = 1:8;
    ax.XTick = 1:epoch_no;
%     Set TickLabels;
    ax.YTickLabel = {'F3','F4','P3','P4','Left','Right','Frontal','Parietal'};
    colormap(my_colormap);
    if sub == 1
        title('Artifact %-age per epoch, (Subjects 1-4, postdrug)');
    elseif sub == all_subs
        xlabel("epoch");
    end
%     colorbar('Ticks',[0.1,0.2,0.5],...
%              'TickLabels',{'10% (wPLI, NC, MFDFA)','20% (ASI)','50% (rest)'});
end
set(gcf, 'Position', get(0, 'Screensize'));
bottom_plot_pos = get(subplot(all_subs,2, all_subs*2),'Position');
second_bottom_plot_pos = get(subplot(all_subs,2, (all_subs-1)*2),'Position');
space_between_plots = second_bottom_plot_pos(2) - (bottom_plot_pos(2) + bottom_plot_pos(4));
colorbar('Position', [bottom_plot_pos(1)+bottom_plot_pos(3)+0.01, ...
                        bottom_plot_pos(2), ...
                        bottom_plot_pos(4)/5, ...
                        (bottom_plot_pos(4)*all_subs + space_between_plots*(all_subs-1))], ...
        'Ticks',[0, 0.01, 0.1,0.2,0.5],...
        'TickLabels',{'', '1% (MFDFA)', '10% (wPLI, NC)','20% (ASI)','50% (rest)'});

savefig(f, fullfile(save_path, 'APT.fig'));
saveas(f, fullfile(save_path, 'APT.png'));
close all;
end

