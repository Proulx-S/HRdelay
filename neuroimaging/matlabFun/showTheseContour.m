function showTheseContour(hLine,ind)

if ~exist('ind','var') || isempty(ind)
    set(hLine,'Visible','on')
else
    set(hLine,'Visible','off')
    set(hLine(ind),'Visible','on')
end
