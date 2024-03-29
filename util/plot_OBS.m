figure
semilogy(necent)
hold on
semilogy(nedist)
semilogy(nedistg)

ax = gca;
grid on
set(get(ax,'xlabel'),'fontsize',14);
set(get(ax,'ylabel'),'fontsize',14);
set(get(ax,'zlabel'),'fontsize',14);
set(get(ax,'title'),'fontsize',14);
set(ax,'fontsize',14);
xlabel('t')
ylabel('estimation error norm')
legend('centralized','distributed','Gramian based distributed','Location','NorthEast')

save('necent.dat','necent','-ascii');
save('nedist.dat','nedist','-ascii');
save('nedistg.dat','nedistg','-ascii');