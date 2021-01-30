clear all
close all

x = load('dp.csv');


plot(smooth(x,25),'-r','linewidth',3)

hold on
for d = [1, 49, 174, 247, 365, 420]
    
    plot([d,d],[0,1],'--k')
    hold on
end

plot([1,420],[0.6, 0.6], '--k','linewidth',2)

set(gca,'fontsize',18)
xlim([1,420])

ylabel('Disorder')
xlabel('Residue')


set(gcf,'PaperUnits','inches','PaperPosition',[0 0 7 3])
print('-painters', '-dpdf', '-r300', 'disorder_plot.pdf');



