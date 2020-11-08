function out = runSVM_patAct4(d,p,verbose,rep,doPerm)

if ~exist('verbose','var')
    verbose = 0;
end
if ~exist('doPerm','var')
    doPerm = 0;
end

%% Permute data if needed
if doPerm
    d.xData = d.xData(randperm(size(d.xData,1)),:);
end

%% Define k-folding
if ischar(p.k)
    switch p.k
        %         case {'auto','autoCmplt'}
        %             k = length(d.crossVal)/2;
        case 'autoRun'
            d.crossVal = defineCrossVal(d,p);
            k = length(d.crossVal)/2/str2double(p.split(1));
        otherwise
            error('');
    end
else
    error('')
    k = p.k;
    p.k = 'randK';
end


%% Loop over cross-validation folds
dOrig = d;
% out.hitRate = nan(k,1);
% out.hitRate_d = nan(k,1);
% out.hitRate_pat = nan(k,1);
% out.hitRate_patRel = nan(k,1);
% 
% out.pat = nan(p.nObs,2);
% out.patRel = nan(p.nObs,2);

out.crossInfo = cell(length(p.trainInfo),length(p.testInfo));
out.crossInfo_dim1Train = p.trainInfo;
out.crossInfo_dim2Test = p.testInfo;
for trainInd = 1:length(p.trainInfo)
    for testInd = 1:length(p.testInfo)
%         display([num2str((trainInd-1)*length(p.testInfo)+testInd) '/' num2str(length(p.trainInfo)*length(p.testInfo))])
        hitRate = nan(k,1);
        decision = nan(p.nObs,1);
        filterdData.f = nan(p.nObs,1);
        filterdData.f1 = nan(p.nObs,1);
        filterdData.f2 = nan(p.nObs,1);
        for fold = 1:k
            [hitRate(fold,1),curDecision,teInd,fData] = trainNtest5(dOrig,p,fold,k,rep,verbose,p.trainInfo{trainInd},p.testInfo{testInd});
            decision(teInd) = curDecision;
            filterdData.f(teInd) = fData.f;
            filterdData.f1(teInd) = fData.f1;
            filterdData.f2(teInd) = fData.f2;
        end
        out.crossInfo{trainInd,testInd}.hitRate = hitRate;
        out.crossInfo{trainInd,testInd}.d = decision;
        out.crossInfo{trainInd,testInd}.f = filterdData.f;
        out.crossInfo{trainInd,testInd}.f1 = filterdData.f1;
        out.crossInfo{trainInd,testInd}.f2 = filterdData.f2;
    end
end

%Recompile
hitRate = nan([size(out.crossInfo) k]);
dVal = nan([size(out.crossInfo) p.nObs]);
f = nan([size(out.crossInfo) p.nObs]);
f1 = nan([size(out.crossInfo) p.nObs]);
f2 = nan([size(out.crossInfo) p.nObs]);
for trainInd = 1:length(p.trainInfo)
    for testInd = 1:length(p.testInfo)
        hitRate(trainInd,testInd,:) = out.crossInfo{trainInd,testInd}.hitRate;
        dVal(trainInd,testInd,:) = out.crossInfo{trainInd,testInd}.d;
        f(trainInd,testInd,:) = out.crossInfo{trainInd,testInd}.f;
        f1(trainInd,testInd,:) = out.crossInfo{trainInd,testInd}.f1;
        f2(trainInd,testInd,:) = out.crossInfo{trainInd,testInd}.f2;
    end
end
out.crossInfo = [];
out.crossInfo.hitRate = hitRate;
out.crossInfo.d = dVal;
out.crossInfo.f = f;
out.crossInfo.f1 = f1;
out.crossInfo.f2 = f2;













function [hitRate,decision_values,teInd,filterdData] = trainNtest5(dOrig,p,fold,k,rep,verbose,trainInfo,testInfo)
%% Train
% Prepare train data
[tr,~,teInd1,teInd2] = prepareDataAtFold(dOrig,p,trainInfo,fold,k,rep,verbose);
trFinal = finalPrep(tr,trainInfo);

teInd = [teInd1; teInd2];

% Train SVM
if ~strcmp(p.kernel,'linear')
    error('did not code for that')
end
% svmStruct = svmtrain(trFinal.data,trFinal.label,'kernel_function',p.kernel,'autoscale',0,'method','SMO','boxconstraint',p.C,'options',statset('MaxIter',1500000,'Display','off'));
svmStruct = svmtrain(trFinal.label,trFinal.data,'-s 0 -t 0 -q');

%% Test
% Prepare test data
[~,te,~,~,~,tePartNorm] = prepareDataAtFold(dOrig,p,testInfo,fold,k,rep,verbose);
teFinal = finalPrep(te,testInfo);
tePartNorm = finalPrep(tePartNorm,'mp');

[~, hitRate, decision_values] = svmpredict4(teFinal.label,teFinal.data,svmStruct,trFinal.data);
hitRate = hitRate/100;

filterdData = wFilter(tePartNorm.data,svmStruct,trFinal.data);




function [tr,te,teInd1,teInd2,trPartNorm,tePartNorm] = prepareDataAtFold(d,p,keepInfo,fold,k,rep,verbose)


if verbose
    tic
    display(['Fold ' num2str(fold) '/' num2str(k) '; Rep ' num2str(rep) '/' num2str(p.repeat)])
end

% Initiate some stuff

%     wFold = w(:,:,:,fold);
%     Afold = A(:,:,:,fold);

% Break complex data down to amp and delay
d2.delay(:,:,1) = angle(d.xData(1:end/2,:));
d2.delay(:,:,2) = angle(d.xData(end/2+1:end,:));
d2.delay(:,:,3) = angle(d.normData);
d2.amp(:,:,1) = abs(d.xData(1:end/2,:));
d2.amp(:,:,2) = abs(d.xData(end/2+1:end,:));
d2.amp(:,:,3) = abs(d.normData);
% switch keepInfo
%     case {'m','p','mp'}
%     case {'r','i'}
%         %remove voxel mean phase estimated from the plaid
%         d2.delay(:,:,1) = wrapToPi(d2.delay(:,:,1)-repmat(circ_mean(d2.delay(:,:,3),[],1),[size(d2.delay,1) 1]));
%         d2.delay(:,:,2) = wrapToPi(d2.delay(:,:,2)-repmat(circ_mean(d2.delay(:,:,3),[],1),[size(d2.delay,1) 1]));
%         d2.delay(:,:,3) = wrapToPi(d2.delay(:,:,3)-repmat(circ_mean(d2.delay(:,:,3),[],1),[size(d2.delay,1) 1]));
%         switch keepInfo
%             case 'r'
%                 %set amplitude to the real value (amplitude at the mean delay)
%                 [d2.amp,~] = pol2cart(d2.delay,d2.amp);
%             case 'i'
%                 %set amplitude to the mag-normalized imaginary value (delay relative to mean delay)
%                 [~,d2.amp] = pol2cart(d2.delay,1);
%         end
%     otherwise
%         error('double-check that')
% end
d2.sessionLabel(:,:,1) = d.sessionLabel(1:end/2,:);
d2.sessionLabel(:,:,2) = d.sessionLabel(end/2+1:end,:);
d2.sessionLabel(:,:,3) = d.sessionLabel(1:end/2,:);
d2.label(:,:,1) = d.label(1:end/2,:);
d2.label(:,:,2) = d.label(end/2+1:end,:);
d2.label(:,:,3) = ones(size(d.label(1:end/2,:))).*3;
if ~strcmp(p.k,'randK')
    d2.crossVal(:,:,1) = d.crossVal(1:end/2,:);
    d2.crossVal(:,:,2) = d.crossVal(end/2+1:end,:);
    d2.crossVal(:,:,3) = nan(size(d.crossVal(1:end/2,:)));
end

switch p.k
    %         case {'auto','autoCmplt'}
    %             error('teInd and trInd badly defined')
    %             teInd = squeeze(d2.crossVal(fold,1,1:2));
    %             d2tr = d2;
    %             allFields = fields(d2tr);
    %             for i = 1:length(allFields)
    %                 d2te.(allFields{i}) = d2tr.(allFields{i})(teInd,:,:);
    %                 d2tr.(allFields{i})(teInd,:,:) = [];
    %             end
    case 'autoRun'
        teInd1 = d2.crossVal(:,1,1)==fold;
        teInd2 = d2.crossVal(:,1,2)==fold;
        allFields = fields(d2);
        for i = 1:length(allFields)
            d2te.(allFields{i})(:,:,1) = d2.(allFields{i})(teInd1,:,1);
            d2te.(allFields{i})(:,:,2) = d2.(allFields{i})(teInd2,:,2);
            tmp1 = d2.(allFields{i})(~teInd1,:,1);
            tmp2 = d2.(allFields{i})(~teInd2,:,2);
            d2tr.(allFields{i}) = cat(3,tmp1,tmp2); clear tmp1 tmp2
        end
        %         case 'sessBal'
        %             error('double-check that')
        %         case 'randK'
        %             % Define train and test data (purely random selection of one sample per session for the testing set, puts more weights on sessions with fewer sample)
        %             sessions = unique(d2.sessionLabel(:,:,1));
        %             for i = 1:length(sessions)
        %                 curInd = find(d2.sessionLabel(:,:,1)==sessions(i));
        %                 teInd(i) = curInd(randperm(length(curInd),1));
        %             end
        %             d2tr = d2;
        %             allFields = fields(d2tr);
        %             for i = 1:length(allFields)
        %                 d2te.(allFields{i}) = d2tr.(allFields{i})(teInd,:,:);
        %                 d2tr.(allFields{i})(teInd,:,:) = [];
        %             end
    otherwise
        error('double-check that')
end


if p.regSession
    error('double-check that')
    %estimate session effect from training set
    [betaHat,X] = estimateSessionEffect(d2tr.amp,d2tr.sessionLabel);
    firsSessInd = size(X,2)-(length(unique(d2tr.sessionLabel))-1)+1;
    sessionEffect = betaHat(:,firsSessInd:end)';
    
    %remove session effect
    for sess = 2:length(unique(d2tr.sessionLabel))
        %from training set
        sessInd = d2tr.sessionLabel(:,:,1)==sess;
        d2tr.amp(sessInd,:,:) = d2tr.amp(sessInd,:,:)-repmat(sessionEffect,[length(find(sessInd)) 1 3]);
        %from testing set
        sessInd = d2te.sessionLabel(:,:,1)==sess;
        d2te.amp(sessInd,:,:) = d2te.amp(sessInd,:,:)-repmat(sessionEffect,[length(find(sessInd)) 1 3]);
    end
end



%% Normalize
allFields = fields(d2tr);
for i = 1:length(allFields)
    tr.(allFields{i}) = cat(1,d2tr.(allFields{i})(:,:,1),d2tr.(allFields{i})(:,:,2));
end
allFields = fields(d2te);
for i = 1:length(allFields)
    te.(allFields{i}) = cat(1,d2te.(allFields{i})(:,:,1),d2te.(allFields{i})(:,:,2));
end

%Subtract mean delay
phase_shift = circ_mean(tr.delay,[],1)+pi/2;
tr.delay = wrapToPi(tr.delay-repmat(phase_shift,size(tr.delay,1),1));
te.delay = wrapToPi(te.delay-repmat(phase_shift,size(te.delay,1),1));
%Scale amp to mean of one
amp_scale = mean(tr.amp,1);
tr.amp = tr.amp./repmat(amp_scale,[size(tr.amp,1) 1]);
te.amp = te.amp./repmat(amp_scale,[size(te.amp,1) 1]);

trPartNorm = tr;
tePartNorm = te;

switch keepInfo
    case 'r'
        %Take real part
        [tr.amp,~] = pol2cart(tr.delay,tr.amp);
        [te.amp,~] = pol2cart(te.delay,te.amp);
        %Set delay to 0
        tr.delay = zeros(size(tr.delay));
        te.delay = zeros(size(te.delay));
%         %Z-score
%         amp_shift = mean(tr.amp,1);
%         tr.amp = tr.amp-repmat(amp_shift,size(tr.amp,1),1);
%         te.amp = te.amp-repmat(amp_shift,size(te.amp,1),1);
%         amp_scale = std(tr.amp,[],1);
%         tr.amp = tr.amp./repmat(amp_scale,size(tr.amp,1),1);
%         te.amp = te.amp./repmat(amp_scale,size(te.amp,1),1);
    case 'i'
        %Take imaginary part
        [~,tr.amp] = pol2cart(tr.delay,tr.amp);
        [~,te.amp] = pol2cart(te.delay,te.amp);
        %Set delay to pi/2
        tr.delay = ones(size(tr.delay))*pi/2;
        te.delay = ones(size(te.delay))*pi/2;
        %         %Z-score
        %         amp_shift = mean(tr.amp,1);
        %         tr.amp = tr.amp-repmat(amp_shift,size(tr.amp,1),1);
        %         te.amp = te.amp-repmat(amp_shift,size(te.amp,1),1);
        %         amp_scale = std(tr.amp,[],1);
        %         tr.amp = tr.amp./repmat(amp_scale,size(tr.amp,1),1);
        %         te.amp = te.amp./repmat(amp_scale,size(te.amp,1),1);
    case 'iN'
        %Remove amp
        tr.amp = ones(size(tr.amp));
        te.amp = ones(size(te.amp));
        %Take imaginary part
        [~,tr.amp] = pol2cart(tr.delay,tr.amp);
        [~,te.amp] = pol2cart(te.delay,te.amp);
        %Set delay to pi/2
        tr.delay = ones(size(tr.delay))*pi/2;
        te.delay = ones(size(te.delay))*pi/2;
        %         %Z-score
        %         amp_shift = mean(tr.amp,1);
        %         tr.amp = tr.amp-repmat(amp_shift,size(tr.amp,1),1);
        %         te.amp = te.amp-repmat(amp_shift,size(te.amp,1),1);
        %         amp_scale = std(tr.amp,[],1);
        %         tr.amp = tr.amp./repmat(amp_scale,size(tr.amp,1),1);
        %         te.amp = te.amp./repmat(amp_scale,size(te.amp,1),1);
    case 'm'
        %Remove delay
        tr.delay = zeros(size(tr.delay));
        te.delay = zeros(size(te.delay));
%         %Z-score)
%         amp_shift = mean(tr.amp,1);
%         tr.amp = tr.amp-repmat(amp_shift,size(tr.amp,1),1);
%         te.amp = te.amp-repmat(amp_shift,size(te.amp,1),1);
%         amp_scale = std(tr.amp,[],1);
%         tr.amp = tr.amp./repmat(amp_scale,size(tr.amp,1),1);
%         te.amp = te.amp./repmat(amp_scale,size(te.amp,1),1);
    case 'p'
        %Remove amp
        tr.amp = ones(size(tr.amp));
        te.amp = ones(size(te.amp));
    case 'mp'
%         %Amp scale only (equivalent to image scaling only)
%         amp_scale = mean(tr.amp,1);
%         tr.amp = tr.amp./repmat(amp_scale,[size(tr.amp,1) 1]);
%         te.amp = te.amp./repmat(amp_scale,[size(te.amp,1) 1]);
%         %Image Shift and Scale
%         [trR,trI] = pol2cart(tr.delay,tr.amp);
%         [teR,teI] = pol2cart(te.delay,te.amp);
%         %Shift
%         %real
%         trR_shift = mean(trR,1);
%         trR = trR-repmat(trR_shift,size(trR,1),1);
%         teR = teR-repmat(trR_shift,size(teR,1),1);
%         %imag
%         trI_shift = mean(trI,1);
%         trI = trI-repmat(trI_shift,size(trI,1),1);
%         teI = teI-repmat(trI_shift,size(teI,1),1);
%         %Scale
%         tr_scale = mean([std(trR,[],1); std(trI,[],1)],1);
%         trR = trR./repmat(tr_scale,size(trR,1),1); trI = trI./repmat(tr_scale,size(trI,1),1);
%         teR = teR./repmat(tr_scale,size(teR,1),1); teI = teI./repmat(tr_scale,size(teI,1),1);
%         
%         [tr.delay,tr.amp] = cart2pol(trR,trI);
%         [te.delay,te.amp] = cart2pol(teR,teI);
%         clear trR trI tr_shift tr_scale teR teI te_shift te_scale
    otherwise
        error('did not code that')
end

if p.doZ
    %Z-score (in the complex domain)
    %train
    [X,Y] = pol2cart(tr.delay,tr.amp);
    
    Xshift = mean(X,1);
    if ~all(Xshift==0); X = X-repmat(Xshift,size(X,1),1); end
    Xscale = std(X,[],1);
    if ~all(Xscale==0); X = X./repmat(Xscale,size(X,1),1); end;
    
    Yshift = mean(Y,1);
    if ~all(Yshift==0); Y = Y-repmat(Yshift,size(Y,1),1); end
    Yscale = std(Y,[],1);
    if ~all(Yscale==0); Y = Y./repmat(Yscale,size(Y,1),1); end;
    
    [tr.delay,tr.amp] = cart2pol(X,Y);
    
    %test
    [X,Y] = pol2cart(te.delay,te.amp);
    
    if ~all(Xshift==0); X = X-repmat(Xshift,size(X,1),1); end
    if ~all(Xscale==0); X = X./repmat(Xscale,size(X,1),1); end;
    
    if ~all(Yshift==0); Y = Y-repmat(Yshift,size(Y,1),1); end
    if ~all(Yscale==0); Y = Y./repmat(Yscale,size(Y,1),1); end;
    
    [te.delay,te.amp] = cart2pol(X,Y);
end





function xData = finalPrep(data,keepInfo)
% switch keepInfo
%     case {'m','r','i','iN'}
%         xData.data = data.amp;
%         xData.label = data.label;
%     case {'p','mp'}
%         [X,Y] = pol2cart(data.delay,data.amp);
%         xData.data = [X Y];
%         xData.label = data.label;
% %     case 'mp'
% %         [X,Y] = pol2cart(data.delay,data.amp);
% %         xData.data = [X Y];
% %         xData.label = data.label;
%     otherwise
%         error('xx')
% end
[X,Y] = pol2cart(data.delay,data.amp);
xData.data = [X Y];
xData.label = data.label;



function [predictedLabel,accuracy,d,out] = svmpredict4(y,x,svmStruct,xTr)
nTr = size(xTr,1);

%% Get w
% w
w = zeros(1,size(svmStruct.SVs,2));
for i = 1:size(svmStruct.SVs,1)
    w = w + svmStruct.sv_coef(i)*svmStruct.SVs(i,:);
end
% w1
ind = find(svmStruct.sv_indices<=nTr/2);
w1 = zeros(1,size(svmStruct.SVs,2));
for i = 1:length(ind)% size(svmStruct.SVs,1)/2
    w1 = w1 + svmStruct.sv_coef(ind(i))*svmStruct.SVs(ind(i),:);
end
% w2
ind = find(svmStruct.sv_indices>nTr/2);
w2 = zeros(1,size(svmStruct.SVs,2));
for i = 1:length(ind)%size(svmStruct.SVs,1)/2+1:size(svmStruct.SVs,1)
    w2 = w2 + svmStruct.sv_coef(ind(i))*svmStruct.SVs(ind(i),:);
end

%% Get a
% a
s = w*xTr';
% a = w*cov(xTr)/cov(s);
a = nan(size(w));
for dim = 1:size(xTr,2)
    tmp = cov(xTr(:,dim),s);
    a(dim) = tmp(1,2);
end

% a1
s1 = w1*xTr(1:end/2,:)';
% a1 = w1*cov(xTr(1:end/2,:))/cov(s1);
a1 = nan(size(w1));
for dim = 1:length(w2)
    tmp = cov(xTr(1:end/2,dim),s1);
    a1(dim) = tmp(1,2);
end

% a2
s2 = w2*xTr(end/2+1:end,:)';
% a2 = w2*cov(xTr(end/2+1:end,:))/cov(s2);
a2 = nan(size(w2));
for dim = 1:length(w2)
    tmp = cov(xTr(end/2+1:end,dim),s2);
    a2(dim) = tmp(1,2);
end

%% Adapt w to test data if needed
if size(x,2)==size(svmStruct.SVs,2)
elseif size(svmStruct.SVs,2)==size(x,2)*2
    w = sqrt(w(1,1:end/2).^2 + w(1,end/2+1:end).^2);
    w1 = sqrt(w1(1,1:end/2).^2 + w1(1,end/2+1:end).^2);
    w2 = sqrt(w2(1,1:end/2).^2 + w2(1,end/2+1:end).^2);
elseif size(svmStruct.SVs,2)==size(x,2)/2
    w = repmat(w,[1 2]);
    w1 = repmat(w1,[1 2]);
    w2 = repmat(w2,[1 2]);
else
    error('xx')
end

%% Get decision value d
% d
d = nan(size(x,1),1);
for j=1:size(x,1)
    d(j,1) = dot(w,x(j,:)) - svmStruct.rho;
end
% d1
d1 = nan(size(x,1),1);
for j=1:size(x,1)
    d1(j,1) = dot(w1,x(j,:)) - svmStruct.rho/2;
end
% d2
d2 = nan(size(x,1),1);
for j=1:size(x,1)
    d2(j,1) = dot(w2,x(j,:)) - svmStruct.rho/2;
end

%% Predict label
predictedLabel = ones(size(d));
predictedLabel(d<0) = 2;

%% Compute accuracy
accuracy = length(find(y==predictedLabel))/length(y)*100;

%% Compile out
out.w = w;
out.w1 = w1;
out.w2 = w2;
out.a = a;
out.a1 = a1;
out.a2 = a2;
out.d = d;
out.d1 = d1;
out.d2 = d2;



function filteredData = wFilter(x,svmStruct,xTr)

nTr = size(xTr,1);

%% Get w
% w
w = zeros(1,size(svmStruct.SVs,2));
for i = 1:size(svmStruct.SVs,1)
    w = w + svmStruct.sv_coef(i)*svmStruct.SVs(i,:);
end
% w1
ind = find(svmStruct.sv_indices<=nTr/2);
w1 = zeros(1,size(svmStruct.SVs,2));
for i = 1:length(ind)% size(svmStruct.SVs,1)/2
    w1 = w1 + svmStruct.sv_coef(ind(i))*svmStruct.SVs(ind(i),:);
end
% w2
ind = find(svmStruct.sv_indices>nTr/2);
w2 = zeros(1,size(svmStruct.SVs,2));
for i = 1:length(ind)%size(svmStruct.SVs,1)/2+1:size(svmStruct.SVs,1)
    w2 = w2 + svmStruct.sv_coef(ind(i))*svmStruct.SVs(ind(i),:);
end


%% Split w and data in real and imaginary parts
if size(w,2)==size(x,2)
    wX = w(1:end/2);
    wY = w(end/2+1:end);
    w1X = w1(1:end/2);
    w1Y = w1(end/2+1:end);
    w2X = w2(1:end/2);
    w2Y = w2(end/2+1:end);
elseif size(w,2)==size(x,2)/2
    wX = w;
    wY = w;
    w1X = w1;
    w1Y = w1;
    w2X = w2;
    w2Y = w2;
else
    error('xx')
end
xX = x(:,1:end/2);
xY = x(:,end/2+1:end);

%% Apply w filter
% d
fX = nan(size(xX,1),1);
for j=1:size(xX,1)
    fX(j,1) = dot(wX,xX(j,:));
end
fY = nan(size(xY,1),1);
for j=1:size(xY,1)
    fY(j,1) = dot(wY,xY(j,:));
end
filteredData.f = complex(fX,fY);
% d1
f1X = nan(size(xX,1),1);
for j=1:size(xX,1)
    f1X(j,1) = dot(w1X,xX(j,:));
end
f1Y = nan(size(xY,1),1);
for j=1:size(xY,1)
    f1Y(j,1) = dot(w1Y,xY(j,:));
end
filteredData.f1 = complex(f1X,f1Y);
% d2
f2X = nan(size(xX,1),1);
for j=1:size(xX,1)
    f2X(j,1) = dot(w2X,xX(j,:));
end
f2Y = nan(size(xY,1),1);
for j=1:size(xY,1)
    f2Y(j,1) = dot(w2Y,xY(j,:));
end
filteredData.f2 = complex(f2X,f2Y);

