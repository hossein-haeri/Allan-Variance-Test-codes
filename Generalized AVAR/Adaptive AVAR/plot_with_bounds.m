function p = plot_with_bounds(x, raw, avg_win, color)
avg = nanmean(reshape([raw'; nan(mod(-numel(raw'),avg_win),1)],avg_win,[]));
std =  nanstd(reshape([raw'; nan(mod(-numel(raw'),avg_win),1)],avg_win,[]));
x = linspace(0,x(end)-x(1),numel(avg));
hold on
patch([x fliplr(x)], [avg-std fliplr(avg+std)], color, 'EdgeColor', 'none', 'FaceAlpha',.2);
p = plot(x', avg, 'Color', color, 'LineWidth', 1);
hold off
end


