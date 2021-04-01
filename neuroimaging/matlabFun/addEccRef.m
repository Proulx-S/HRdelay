function addEccRef(d,p)
theta = linspace(-pi,pi,100);
eccRef = d.voxProp.eccTrans(p.featSel.fov.threshVal);
for i = 1:length(eccRef)
    [x,y] = pol2cart(theta,eccRef(i));
%     plot(x,y,'--','Color',[1 1 1].*0.4)
    plot(x,y,'--k')
end
