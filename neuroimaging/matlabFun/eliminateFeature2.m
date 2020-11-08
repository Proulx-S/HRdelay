function ind2keep = eliminateFeature2(tr,trLabel,p,nFeatIn)
global err
autoscaleFlag = false;

%% Loop over the scoring folds
wL = nan(p.RFE_L,size(tr,2));
for Lfold = 1:p.RFE_L
    LsplitInd = [randperm(size(tr,1)/2)'>size(tr,1)/2*(5/20); randperm(size(tr,1)/2)'>size(tr,1)/2*(5/20)];
    
    tr_L = tr(LsplitInd,:);
    trLabel_L = trLabel(LsplitInd);
    
    %% Train the model
    switch p.kernel
        case 'polynomial'
            svmParams = statset('MaxIter',1500000,'Display','off');
            %         try
            svmStruct = svmtrain(tr,trLabel,'kernel_function',p.kernel,'polyorder',p.polyorder,'autoscale',autoscaleFlag,'method','SMO','boxconstraint',p.C,'options',svmParams);
            %         catch err
            %             return
            %         end
        otherwise
            svmParams = statset('MaxIter',1500000,'Display','off');
            %         try
            svmStruct = svmtrain(tr,trLabel,'kernel_function',p.kernel,'autoscale',autoscaleFlag,'method','SMO','boxconstraint',p.C,'options',svmParams);
            %         catch err
            %             warning(err.identifier,err.message)
            %             return
            %         end
    end
    
    %% Get filter weigths
    alpha = zeros(size(tr_L,1),1);
    alpha(svmStruct.SupportVectorIndices,1) = svmStruct.Alpha;
    y = (trLabel_L-1)*2-1;
    x = tr_L;
    w = zeros(1,size(x,2));
    for i = 1:size(x,1)
        w = w + alpha(i)*y(i)*x(i,:);
    end
    
    %% Transform the weights to scallar when necessary
    w = abs(w);
%     switch p.learningVariable % 'r' 'i' 'ri' 'm' 'p' 'mp'
%         case 'ri'
%             w = sqrt(w(1:end/2).^2+w(end/2+1:end).^2);
%         case {'r','i'}
%             w = abs(w);
%         case {'m','p'}
%             w = abs(w);
%         case {'mp'}
%             error('not implemented')
%         otherwise
%             error('p.curParam badly specifie')
%     end
    wL(Lfold,:) = w;
end
%% Sum filter weights
w = sum(wL,1)/p.RFE_L;
[~,b] = sort(w,'descend');
ind2keep = b(1:nFeatIn);
