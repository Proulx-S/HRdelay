function rad = sec2rad(sec,period)
rad = nan(size(sec));
for i = 1:size(sec,1)
    for ii = 1:size(sec,2)
        for iii = 1:size(sec,3)
            rad(i,ii,iii) = (sec(i,ii,iii)-period/2)/(period/2)*pi;
        end
    end
end




