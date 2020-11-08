function sec = rad2sec(rad,period)
sec = nan(size(rad));
for i = 1:size(sec,1)
    for ii = 1:size(sec,2)
        for iii = 1:size(sec,3)
            sec(i,ii,iii) = rad(i,ii,iii)/pi*(period/2)+(period/2);
        end
    end
end