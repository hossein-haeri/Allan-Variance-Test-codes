function plot_results(data_t,data_y,t_true,y_true,avar,mse,tau,experiment_name)
%% visualization
set(groot,'defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultAxesFontSize',16);

[min_avar,idx_avar] = min(mean(avar,2));
[min_mse, idx_mse] = min(mean(mse,2));
num_est_samples = 500;
[t_hat, y_hat] = estimate_at(data_t,data_y,tau(idx_avar),num_est_samples);
[t_hat_mse, y_hat_mse] = estimate_at(data_t,data_y,tau(idx_mse),num_est_samples);


c1 = [0.6350 0.0780 0.1840];
c2 = [0.3010 0.7450 0.9330];

fig1 = figure('Position', [10 10 550 350]);
    hold on
%     for i=1:numel(data_t)
%         xline(data_t(i),'LineStyle','-','LineWidth',.5,'Color',[.8 .8 .8]);
%     end
    scatter(data_t,data_y,40,'.','MarkerEdgeAlpha',0.4,'MarkerEdgeColor',[.4 .4 .4]);
    p_true = plot(t_true,y_true,          'LineWidth',1.0,'Color',[.4 .4 .4],'LineStyle','-');
%     p_mse = plot(t_hat_mse,y_hat_mse,'-^','LineWidth',1,'Color',[0.8500, 0.3250, 0.0980]);
%     p_hat = plot(t_hat,y_hat,        '-*','LineWidth',1,'Color',[0, 0.4470, 0.7410]);
    p_mse = plot(t_hat_mse,y_hat_mse,'LineWidth',3,'Color',[c2 0.7]);
    p_hat = plot(t_hat,y_hat,'LineWidth',1.1,'Color',c1);
    p_hat.MarkerSize = 9;
    p_mse.MarkerSize = 9;
    p_hat.MarkerIndices = 1:25:length(y_hat);
    p_mse.MarkerIndices = 1:25:length(y_hat_mse);

%     xline(t_range/2-tau(idx)/2,'LineStyle','--','LineWidth',2,'Color',[.3 .3 .3]);
%     xline(t_range/2+tau(idx)/2,'LineStyle','--','LineWidth',2,'Color',[.3 .3 .3]);
    xlabel('$t [s]$')
    ylabel('$\theta$')
    y_min=min(data_y);
    y_max=max(data_y);
    ylim([y_min y_max+(y_max-y_min)/3.5])
    grid on
    legend('Data','Actual','MSE Est.','AVAR Est.','Location','north','Orientation','horizontal','FontSize',10.5,'EdgeColor',[.4 .4 .4])

fig2 = figure('Position', [400 10 550 350]);
    subplot(2,1,1)
        hold on
        avg_avar = mean(avar,2)';
        std_avar = std(avar,0,2)';
        patch([tau fliplr(tau)], [avg_avar-std_avar fliplr(avg_avar+std_avar)], [.2 .2 .2], 'EdgeColor', 'none', 'FaceAlpha',.1);
        plot(tau, avg_avar,'LineWidth',1.5, 'Color',[.3 .3 .3]);
        scatter(tau(idx_avar),min_avar,60,'s','filled','MarkerEdgeColor',c1,'MarkerFaceColor',c1)
        text(tau(idx_avar)*6/10,6/3*min_avar,['$\tau = ' num2str(round(tau(idx_avar),4)) '$ s'],'Interpreter','latex','FontSize',11)
        xlabel('Window length $\tau [s]$')
%         ylabel('AVAR $\sigma^2_A$','Position',[1.3,0.031315763324153,-1])
        ylabel('AVAR $\sigma^2_A$')
        set(gca,'xscale','log')
        set(gca,'yscale','log')
        grid on
        xlim([tau(1) tau(end)]);
%         title('(a)')
        

    subplot(2,1,2)
        hold on
        avg_mse = mean(mse,2)';
        std_mse = std(mse,0,2)';
        patch([tau fliplr(tau)], [avg_mse-std_mse fliplr(avg_mse+std_mse)], [.2 .2 .2], 'EdgeColor', 'none', 'FaceAlpha',.1);
        plot(tau, avg_mse,'LineWidth',1.5, 'Color',[.3 .3 .3]);
        scatter(tau(idx_mse),min_mse,60,'s','filled','MarkerEdgeColor',c2,'MarkerFaceColor',c2)
        text(round(tau(idx_mse),2)*6/10,6/3*min_mse,['$\tau = ' num2str(round(tau(idx_mse),4)) '$ s'],'Interpreter','latex','FontSize',11)
        xlabel('Window length $\tau [s]$')
%         ylabel('Estimation MSE','Position',[1.3,0.175180343921021])
        ylabel('Estimation MSE')
        set(gca,'xscale','log')
        set(gca,'yscale','log')
        grid on
        xlim([tau(1) tau(end)]);
%         title('(b)')
        
align_Ylabels(fig2)


saveas(fig1,['example_' experiment_name '_data'],'svg')  
saveas(fig2,['example_' experiment_name '_avar'],'svg')   
end







