function plotDensity(x)

x = [real(x) imag(x)];
[~,density,X,Y]=kde2d(x);
figure('WindowStyle','docked');
surf(X,Y,density,'LineStyle','none'), view([0,70])
colormap hot, hold on, alpha(.8)
set(gca, 'color', 'b');

xDensity = interp2(X,Y,density,x(:,1),x(:,2));
scatter3(x(:,1),x(:,2),xDensity,eps,'w.')
ax = gca;
ax.DataAspectRatio = ax.DataAspectRatio([1 1 3]);
ax.PlotBoxAspectRatio = ax.PlotBoxAspectRatio([1 1 3]);
plot3(xlim,[0 0],[1 1].*max(xDensity)/2,'w')
plot3([0 0],ylim,[1 1].*max(xDensity)/2,'w')