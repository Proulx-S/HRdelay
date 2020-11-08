function mi = computeMI(tr,trLabel1,trLabel2)
if isreal(tr)
    isReal = 1;
    % real only
    xData1a = tr(trLabel1,:);
    xData2a = tr(trLabel2,:);
else
    isReal = 0;
    % real+imag
    xData1a = real(tr(trLabel1,:));
    xData2a = real(tr(trLabel2,:));
    xData1b = imag(tr(trLabel1,:));
    xData2b = imag(tr(trLabel2,:));
    if any(all(xData1a==0,1)) || any(all(xData2a==0,1))
        error('no variance on the real part')
    end
    if any(all(xData1b==0,1)) || any(all(xData2b==0,1))
        error('no variance on the imag part')
    end
end


yData = trLabel2;
featSelVala = nan(1,size(xData1a,2));
for vox = 1:size(xData1a,2)
    xData = [xData1a(:,vox); xData2a(:,vox)];
    featSelVala(1,vox) = information(xData',yData');
end
if ~isReal
    featSelValb = nan(1,size(xData1a,2));
    for vox = 1:size(xData1a,2)
        xData = [xData1b(:,vox); xData2b(:,vox)];
        featSelValb(1,vox) = information(xData',yData');
    end
    mi = (featSelVala+featSelValb)./2;
else
    mi = featSelVala;
end
