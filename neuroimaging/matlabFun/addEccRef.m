function [x,y] = addEccRef(voxProp,p)
theta = linspace(-pi,pi,100);
eccRef = voxProp.eccTrans(p.featSel.fov.threshVal);
x = nan(size(theta,2),length(eccRef));
y = nan(size(theta,2),length(eccRef));
for i = 1:length(eccRef)
    [x(:,i),y(:,i)] = pol2cart(theta,eccRef(i));
%     plot(x,y,'--','Color',[1 1 1].*0.4)
    plot(x(:,i),y(:,i),'--k')
end
