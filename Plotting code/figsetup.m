function figsetup(fig_no)
    
    % Format figure settings 
    
    figure(fig_no); hold on;
    set(gcf,'Position',[360 278 500 500])
    ax1 = gca;
    ax1.FontSize = 20;
    ax1.TitleFontSizeMultiplier = 1;
    ax1.LabelFontSizeMultiplier = 1;
    ax1.LineWidth = 2;
    axis square
    ax1.Position = [0.2 0.2 0.6 0.6];
    
    set(gca, 'Layer', 'top')
end