function sPlotLogisticRegressionWithData(x_vals,epts,b_vals,stats,name)

x_axis = (0:max(x_vals));
[rpb,rplo,rphi] = glmval(b_vals, x_axis,'logit',stats);

figure('Name',['Heart ',name,' vs. RP'],'NumberTitle','off');
hold on;

plot(x_axis,rpb,'b','LineWidth',2); % responding function
plot(x_axis,rpb-rplo,'b','LineWidth',1); % low CI curve
plot(x_axis,rpb+rphi,'b','LineWidth',1); % high CI curve
%set(gca,'xminortick','on','yminortick','on');
%set(gca,'box','on');
%xlabel('gEUD'); ylabel('RP probability');

% Load patient values
flg = ~epts;
[medians,means,stdmeans,binlow,binhigh,numcomp,numtotal,betainv84,betainv16]...
    = EventObserved(flg,x_vals,4);
prob = numcomp./numtotal;
% plot
errorbar(medians,prob,max(0,prob-betainv16),max(0,betainv84-prob),'r*','LineWidth',1);
% horizontal error bars

errorbar_x(means,prob,(means-binlow),(binhigh-means),'r');


%ploterr(means,prob,stdmeans,zeros(size(prob)),'r');


ylims = ylim;
ymax = ylims(2);
ylim([0 ymax]);

xlims = xlim;
xmax = xlims(2);
xlim([0 xmax]);

set(gca,'xminortick','on','yminortick','on');
set(gca,'box','on');
xlabel(name); ylabel('Observed rate of RP3');

stat_text = 'B = [%6.3g, %6.3g]\n\np = [%6.3g, %6.3g]';
Pvals = stats.p;
%text(10,0.45,...
text(4,0.45,...
    sprintf(stat_text,[b_vals(1) b_vals(2) Pvals(1) Pvals(2)]),...
    'FontSize',10,'BackgroundColor','w','EdgeColor','k');