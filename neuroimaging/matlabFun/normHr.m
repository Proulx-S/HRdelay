function hr = normHr(hr,hrSin)
interpMethod = 'cubic'; % cubic convolution as implemented in Matlab2020b

% Seperate group mean from between subject variations
hrSin = mean(hrSin(:,:,:),3);
rhoSubj = abs(hrSin);
thetaSubj = angle(hrSin);
rhoGroup = abs(mean(hrSin,2));
thetaGroup = angle(mean(hrSin,2));
theta = -( thetaGroup-thetaSubj );
rho = rhoSubj./rhoGroup;

% Define time
tPts = 12;
deltaSec = 1;
tSec = linspace(0,(tPts-1)*deltaSec,tPts);
deltaSec2 = 0.001;
tPts2 = tPts*(deltaSec/deltaSec2);
tSec2 = linspace(0,(tPts2-1)*deltaSec2,tPts2);

deltaRad2 = deltaSec2/(tPts2*deltaSec2)*2*pi;

ind = 1:size(tSec,2);
indPreppend = size(tSec,2)/2:size(tSec,2);
indAppend = 1:size(tSec,2)/2;
tSecPreppend = -length(indPreppend)*deltaSec:deltaSec:tSec(1)-deltaSec;
tSecAppend = tSec(end)+deltaSec:deltaSec:tSec(end)+length(indAppend)*deltaSec;
indPreppend2 = size(tSec2,2)/2:size(tSec2,2);
indAppend2 = 1:size(tSec2,2)/2;
tSecPreppend2 = -length(indPreppend2)*deltaSec2:deltaSec2:tSec2(1)-deltaSec2;
tSecAppend2 = tSec2(end)+deltaSec2:deltaSec2:tSec2(end)+length(indAppend2)*deltaSec2;

for bRepInd = 1:size(hr,2)
    for wRepInd = 1:numel(hr(1,1,:))
        % Upsample
        tmp = interp1([tSecPreppend tSec tSecAppend],hr([indPreppend ind indAppend],bRepInd,wRepInd),[tSecPreppend2 tSec2 tSecAppend2],interpMethod);
        % Rotate
        tmp = circshift(tmp,round(theta(1,bRepInd)./deltaRad2));
        % Scale
        tmp = tmp./rho(1,bRepInd);
        % Downsample
        tmp = interp1([tSecPreppend2 tSec2 tSecAppend2],tmp,[tSecPreppend tSec tSecAppend],interpMethod);
        tmp = tmp(length(tSecPreppend)+1:length([tSecPreppend tSec]));
        hr(:,bRepInd,wRepInd) = tmp;
    end
end
