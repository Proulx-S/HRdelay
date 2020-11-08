function w = crossInfoAdapt(w,info,p)

if ~strcmp(info{1},info{2})
    if strcmp(info{1},'r') && strcmp(info{2},'mp')
    elseif strcmp(info{1},'r') && strcmp(info{2},'m')
        
    elseif strcmp(info{1},'r') && strcmp(info{2},'p')
        switch p.crossInfoMethod
            case 'rot'
                [wDelay,wAmp] = cart2pol(w(1:end/2),w(end/2+1:end));
                [w(1:end/2),w(end/2+1:end)] = pol2cart(wDelay+pi/2,wAmp);
            case 'proj'
                w = cat(2,zeros(size(w(end/2+1:end))),w(1:end/2));
        end
    elseif strcmp(info{1},'r') && strcmp(info{2},'i')
        switch p.crossInfoMethod
            case 'rot'
                [wDelay,wAmp] = cart2pol(w(1:end/2),w(end/2+1:end));
                [w(1:end/2),w(end/2+1:end)] = pol2cart(wDelay+pi/2,wAmp);
            case 'proj'
                w = cat(2,zeros(size(w(end/2+1:end))),w(1:end/2));
        end
    elseif strcmp(info{1},'mp') && strcmp(info{2},'r')
        w = [w(1:end/2) zeros(size(w(end/2+1:end)))];
    elseif strcmp(info{1},'mp') && strcmp(info{2},'m')        
        w = [w(1:end/2) zeros(size(w(end/2+1:end)))];
    elseif strcmp(info{1},'mp') && strcmp(info{2},'p')
        w = [zeros(size(w(1:end/2))) w(end/2+1:end)];
    elseif strcmp(info{1},'mp') && strcmp(info{2},'i')
        w = [zeros(size(w(1:end/2))) w(end/2+1:end)];
    elseif strcmp(info{1},'m') && strcmp(info{2},'r')
        
    elseif strcmp(info{1},'m') && strcmp(info{2},'mp')
        
    elseif strcmp(info{1},'m') && strcmp(info{2},'p')
        switch p.crossInfoMethod
            case 'rot'
                [wDelay,wAmp] = cart2pol(w(1:end/2),w(end/2+1:end));
                [w(1:end/2),w(end/2+1:end)] = pol2cart(wDelay+pi/2,wAmp);
            case 'proj'
                w = cat(2,zeros(size(w(end/2+1:end))),w(1:end/2));
        end
    elseif strcmp(info{1},'m') && strcmp(info{2},'i')
        switch p.crossInfoMethod
            case 'rot'
                [wDelay,wAmp] = cart2pol(w(1:end/2),w(end/2+1:end));
                [w(1:end/2),w(end/2+1:end)] = pol2cart(wDelay+pi/2,wAmp);
            case 'proj'
                w = cat(2,zeros(size(w(end/2+1:end))),w(1:end/2));
        end
    elseif strcmp(info{1},'p') && strcmp(info{2},'r')
        switch p.crossInfoMethod
            case 'rot'
                [wDelay,wAmp] = cart2pol(w(1:end/2),w(end/2+1:end));
                [w(1:end/2),w(end/2+1:end)] = pol2cart(wDelay-pi/2,wAmp);
            case 'proj'
                w = cat(2,w(end/2+1:end),zeros(size(w(1:end/2))));
        end
    elseif strcmp(info{1},'p') && strcmp(info{2},'mp')
        
    elseif strcmp(info{1},'p') && strcmp(info{2},'m')
        switch p.crossInfoMethod
            case 'rot'
                [wDelay,wAmp] = cart2pol(w(1:end/2),w(end/2+1:end));
                [w(1:end/2),w(end/2+1:end)] = pol2cart(wDelay-pi/2,wAmp);
            case 'proj'
                w = cat(2,w(end/2+1:end),zeros(size(w(1:end/2))));
        end
    elseif strcmp(info{1},'p') && strcmp(info{2},'i')
        
    elseif strcmp(info{1},'i') && strcmp(info{2},'mp')
        
    elseif strcmp(info{1},'i') && strcmp(info{2},'m')
        switch p.crossInfoMethod
            case 'rot'
                [wDelay,wAmp] = cart2pol(w(1:end/2),w(end/2+1:end));
                [w(1:end/2),w(end/2+1:end)] = pol2cart(wDelay-pi/2,wAmp);
            case 'proj'
                w = cat(2,w(end/2+1:end),zeros(size(w(1:end/2))));
        end
    elseif strcmp(info{1},'i') && strcmp(info{2},'p')
        
    elseif strcmp(info{1},'i') && strcmp(info{2},'r')
        switch p.crossInfoMethod
            case 'rot'
                [wDelay,wAmp] = cart2pol(w(1:end/2),w(end/2+1:end));
                [w(1:end/2),w(end/2+1:end)] = pol2cart(wDelay-pi/2,wAmp);
            case 'proj'
                w = cat(2,w(end/2+1:end),zeros(size(w(1:end/2))));
        end
    else
        error('')
    end
end