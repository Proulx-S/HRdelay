function indIn = contour2mask(hLine,ind,hScat)

if ~exist('ind','var') || isempty(ind)
    set(hLine,'Visible','on')
else
    set(hLine,'Visible','off')
    set(hLine(abs(ind)),'Visible','on')
end
if exist('hScat','var')
    indIn = false(size(hScat.XData,2),length(ind));
    for contInd = 1:length(ind)
        indIn(:,contInd) = inpolygon(hScat.XData,hScat.YData,hLine(abs(ind(contInd))).XData,hLine(abs(ind(contInd))).YData);
        if ind(contInd)<0
            indIn(:,contInd) = ~indIn(:,contInd);
        end
    end
    indIn = all(indIn,2);
%     hScat2 = copyobj(hScat,hScat.Parent);
%     hScat2.XData(:,~indIn) = [];
%     hScat2.YData(:,~indIn) = [];
%     hScat.Visible = 'off';
else
    indIn = [];
end