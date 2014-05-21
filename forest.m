function forest(param,texts, mean,low,high,size2,...
    q, q_pval, t2, i2)

subplot(1,3,1)
axis([0,1.5,0,length(texts)+1])
set(gca,'visible','off')

for n=1:length(texts)
    text(-0.5,n,texts{n},'FontSize',18,'interpreter','latex');
    numbers = sprintf('$%6.2f~[%6.2f; %6.2f]$',mean(n),low(n),high(n));
    text(0.36,n,numbers,'FontSize',18,'interpreter','latex');
end
text(0.32,n+1,['\underline{\textbf{\textit{',param,'}}~[95\% CIs]}'],...
    'FontSize',20,'interpreter','latex');




subplot(1,3,[2,3])


plot([low,high]',[1;1]*(1:length(mean)),'k');
axis([min([min(low)-0.1*min(low),-0.1]),max(high)+0.1*max(high),0,length(mean)+1])
hold on
for i=1:length(texts)
plot(mean(i),i,'ko','markerSize',size2(i)+1,'MarkerEdgeColor','k','MarkerFaceColor','k')
end

text(0.55,0.35,['$Q = ',num2str(round(q*100)/100,3),',~p = ',num2str(round(q_pval*1000)/1000,4),'$'],...
    'Unit','normalize',...
    'FontSize',16,'interpreter','latex'); 

text(0.55,0.25,['$\tau^2 = ',num2str(round(t2*100)/100,3),',~I^2 = ',num2str(round(i2*100)/100,3),'$'],...
    'Unit','normalize',...
    'FontSize',16,'interpreter','latex'); 

hold off;
% grid on;

set(gca,'FontSize',20);

set(gca,'YTickLabel',[]);
xlabel(['LKB ',param],'FontSize',22,'interpreter','latex');

ax = axes('Position',[0 0 1 1],'Unit','normalize',...
    'parent',gcf);
plot([0.06 0.905],[0.47 0.47],'k--','LineWidth',1.5);
set(ax,'xlim',[0 1]);
set(ax,'ylim',[0 1]);

set(ax,'visible','off');
end