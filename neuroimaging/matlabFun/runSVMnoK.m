function svm = runSVMnoK(d,p,verbose,ppn,doPerm)
if ~exist('verbose','var')
    verbose = 1;
end
if ~exist('ppn','var')
    ppn = 0;
end
if ~exist('doPerm','var')
    doPerm = 0;
    p.doPerm = 0;
else
    p.doPerm = doPerm;
end

switch p.algorithm
    case 'runSVM'
        hitRate = nan(p.repeat,p.nFeatLev,p.nFeatLev);
    case 'runSVM_RFE'
        hitRate = nan(p.repeat,p.nFeatLev);
    case 'runSVM_patAct'
        if ischar(p.k)
            switch p.k
                %         case {'auto','autoCmplt'}
                %             k = length(d.crossVal)/2;
                case 'autoRun'
                    tmp = defineCrossVal(d,p);
                    k = length(tmp)/2/str2double(p.split(1));
                otherwise
                    error('');
            end
        else
            error('')
            k = p.k;
            p.k = 'randK';
        end
        
        out.hitRate = nan(p.repeat,k);
        out.hitRate_d = nan(p.repeat,k);
        out.hitRate_pat = nan(p.repeat,k);
        out.hitRate_patRel = nan(p.repeat,k);
        out.pat = nan(p.repeat,p.nObs,2);
        out.patRel = nan(p.repeat,p.nObs,2);
end


if ppn
    %% Run SVM with parallel matlab
    numworkers = ppn-1;
    if verbose
        display(['Running SVM at ' num2str(length(p.featLevList)) '*' num2str(length(p.featLevList)) '/2 feature selection levels'])
        display(['for ' num2str(p.repeat) ' repetitions on ' num2str(numworkers) ' workers'])
    end
    
    
    clust = parcluster('local');
    if isempty(getenv('OS'))
        tmpJobStorageLocation = fullfile(clust.JobStorageLocation,['clust_' datestr(now,'yyyymmddHHMMSSFFF') num2str(round(rand(1,1)*10)) num2str(round(rand(1,1)*10)) num2str(round(rand(1,1)*10))]);
        mkdir(tmpJobStorageLocation)
        clust.JobStorageLocation = tmpJobStorageLocation;
    end
    
    poolobj = gcp('nocreate'); % If no pool, do not create new one.
    if isempty(poolobj)
        poolobj = parpool(clust,numworkers);
    elseif poolobj.NumWorkers~=numworkers
        delete(poolobj);
        poolobj = parpool(clust,numworkers);
    else
        display(['parpool already created with ' num2str(numworkers) ' workers'])
    end
    
    % Run SVM
    switch p.algorithm
%         case 'runSVM'
%             parfor rep = 1:p.repeat
%                 hitRate(rep,:,:) = runSVM6(d,p,verbose,rep,doPerm);
%             end
%         case 'runSVM_RFE'
%             parfor rep = 1:p.repeat
%                 hitRate(rep,:) = runSVM6(d,p,verbose,rep,doPerm);
%             end
%         case 'runSVM_patAct2'
%             error('double-check that')
%             parfor rep = 1:p.repeat
%                 outRep{rep} = runSVM_patAct3(d,p,0,rep,doPerm);
%             end
%             for rep = 1:p.repeat
%                 %Compile
%                 out.crossInfo.hitRate(rep,:,:) = outRep{rep}.crossInfo.hitRate;
%                 out.crossInfo.d(rep,:,:,:) = outRep{rep}.crossInfo.d;
%                 out.crossInfo_dim1Train = outRep{rep}.crossInfo_dim1Train;
%                 out.crossInfo_dim2Test = outRep{rep}.crossInfo_dim2Test;
%             end
%             clear outRep
%         case 'runSVM_patAct4'
%             parfor rep = 1:p.repeat
%                 outRep{rep} = runSVM_patAct4(d,p,0,rep,doPerm);
%             end
%             for rep = 1:p.repeat
%                 %Compile
%                 out.crossInfo.hitRate(rep,:,:,:) = outRep{rep}.crossInfo.hitRate;
%                 out.crossInfo.d(rep,:,:,:) = outRep{rep}.crossInfo.d;
%                 out.crossInfo.f(rep,:,:,:) = outRep{rep}.crossInfo.f;
%                 out.crossInfo.f1(rep,:,:,:) = outRep{rep}.crossInfo.f1;
%                 out.crossInfo.f2(rep,:,:,:) = outRep{rep}.crossInfo.f2;
%                 out.crossInfo_dim1Train = outRep{rep}.crossInfo_dim1Train;
%                 out.crossInfo_dim2Test = outRep{rep}.crossInfo_dim2Test;
%             end
%             clear outRep
%         case 'runSVM_patAct5'
%             parfor rep = 1:p.repeat
%                 outRep{rep} = runSVM_patAct5(d,p,0,rep,doPerm);
%             end
%             for rep = 1:p.repeat
%                 %Compile
%                 out.crossInfo.hitRate(rep,:,:,:) = outRep{rep}.crossInfo.hitRate;
%                 out.crossInfo.d(rep,:,:,:) = outRep{rep}.crossInfo.d;
%                 out.crossInfo.f(rep,:,:,:) = outRep{rep}.crossInfo.f;
%                 out.crossInfo.f1(rep,:,:,:) = outRep{rep}.crossInfo.f1;
%                 out.crossInfo.f2(rep,:,:,:) = outRep{rep}.crossInfo.f2;
%                 out.crossInfo_dim1Train = outRep{rep}.crossInfo_dim1Train;
%                 out.crossInfo_dim2Test = outRep{rep}.crossInfo_dim2Test;
%             end
%             clear outRep
        case 'runSVM_patAct6'
            parfor rep = 1:p.repeat
                outRep{rep} = runSVM_patAct6_noK(d,p,0,rep);;
            end
            for rep = 1:p.repeat
                %Compile
%                 out.crossInfo.hitRate(rep,:,:,:) = outRep{rep}.crossInfo.hitRate;
%                 out.crossInfo.d(rep,:,:,:) = outRep{rep}.crossInfo.d;
                out.crossInfo.w(rep,:,:,:) = outRep{rep}.w;
                out.crossInfo.a(rep,:,:,:) = outRep{rep}.a;
%                 out.crossInfo.f(rep,:,:,:) = outRep{rep}.crossInfo.f;
%                 out.crossInfo.f1(rep,:,:,:) = outRep{rep}.crossInfo.f1;
%                 out.crossInfo.f2(rep,:,:,:) = outRep{rep}.crossInfo.f2;
                out.crossInfo_dim1Train = p.trainInfo;
%                 out.crossInfo_dim2Test = outRep{rep}.crossInfo_dim2Test;
%                 switch p.infoComb
%                     case 'layeredSVM'
%                         out.newRes.hitRate(rep,:) = outRep{rep}.newRes.hitRate;
%                         out.newRes.d(rep,:,:) = outRep{rep}.newRes.d;
%                         out.newRes.svmBias(rep,:,:) = outRep{rep}.newRes.svmBias;
%                         out.newRes.w(rep,:,:) = squeeze(outRep{rep}.newRes.w);
%                         out.newRes.a(rep,:,:) = squeeze(outRep{rep}.newRes.a);
%                         out.newRes.info(rep,:) = outRep{rep}.newRes.info;
%                     case 'layeredSVMwCross'
%                         out.combInfo.hitRate(rep,:) = outRep{rep}.combInfo.hitRate;
%                         out.combInfo.d(rep,:,:) = outRep{rep}.combInfo.decision';
%                     otherwise
%                         error('missing something');
%                 end
            end
            clear outRep
    end
    if isempty(getenv('OS'))
        delete(poolobj)
        rmdir(tmpJobStorageLocation,'s');
    end
else
    %% Run SVM with regular matlab
    for rep = 1:p.repeat
        if verbose
            time1 = tic;
            display(['Starting repetition ' num2str(rep) '/' num2str(p.repeat)])
        end
        
        % Run SVM
        switch p.algorithm
%             case 'runSVM'
%                 hitRate(rep,:,:) = runSVM6(d,p,verbose,rep,doPerm);
%             case 'runSVM_RFE'
%                 keyboard
%                 hitRate(rep,:) = runSVM6(d,p,verbose,rep,doPerm);
%             case 'runSVM_patAct'
%                 outRep = runSVM_patAct2(d,p,verbose,rep,doPerm);
%                 %Compile
%                 allFields = fields(outRep);
%                 for i = 1:length(allFields)
%                     if strfind(allFields{i},'hitRate')
%                         out.(allFields{i})(rep,:) = outRep.(allFields{i});
%                     elseif strfind(allFields{i},'pat')
%                         out.(allFields{i})(rep,:,:) = outRep.(allFields{i});
%                     else
%                         error('something wrong here')
%                     end
%                 end
%                 clear outRep
%             case 'runSVM_patAct2'
%                 error('double-check that')
%                 outRep = runSVM_patAct3(d,p,0,rep,doPerm);
%                 
%                 %Compile
%                 out.crossInfo.hitRate(rep,:,:) = outRep.crossInfo.hitRate;
%                 out.crossInfo.d(rep,:,:,:) = outRep.crossInfo.d;
%                 out.crossInfo_dim1Train = outRep.crossInfo_dim1Train;
%                 out.crossInfo_dim2Test = outRep.crossInfo_dim2Test;
%                 clear outRep
%             case 'runSVM_patAct4'
%                 outRep = runSVM_patAct4(d,p,0,rep,doPerm);
%                 %Compile
%                 out.crossInfo.hitRate(rep,:,:,:) = outRep.crossInfo.hitRate;
%                 out.crossInfo.d(rep,:,:,:) = outRep.crossInfo.d;
%                 out.crossInfo.f(rep,:,:,:) = outRep.crossInfo.f;
%                 out.crossInfo.f1(rep,:,:,:) = outRep.crossInfo.f1;
%                 out.crossInfo.f2(rep,:,:,:) = outRep.crossInfo.f2;
%                 out.crossInfo_dim1Train = outRep.crossInfo_dim1Train;
%                 out.crossInfo_dim2Test = outRep.crossInfo_dim2Test;
%                 clear outRep
%             case 'runSVM_patAct5'
%                 outRep = runSVM_patAct5(d,p,0,rep,doPerm);
%                 %Compile
%                 out.crossInfo.hitRate(rep,:,:,:) = outRep.crossInfo.hitRate;
%                 out.crossInfo.d(rep,:,:,:) = outRep.crossInfo.d;
%                 out.crossInfo.f(rep,:,:,:) = outRep.crossInfo.f;
%                 out.crossInfo.f1(rep,:,:,:) = outRep.crossInfo.f1;
%                 out.crossInfo.f2(rep,:,:,:) = outRep.crossInfo.f2;
%                 out.crossInfo_dim1Train = outRep.crossInfo_dim1Train;
%                 out.crossInfo_dim2Test = outRep.crossInfo_dim2Test;
%                 clear outRep
            case 'runSVM_patAct6'
%                 p.trainInfo = {'p'}; p.testInfo = {'r'};
                outRep = runSVM_patAct6_noK(d,p,0,rep);
                %Compile
%                 out.crossInfo.hitRate(rep,:,:,:) = outRep.crossInfo.hitRate;
%                 out.crossInfo.d(rep,:,:,:) = outRep.crossInfo.d;
                out.crossInfo.w(rep,:,:) = outRep.w; % rep x info x vox
                out.crossInfo.a(rep,:,:) = outRep.a; % rep x info x vox
%                 out.crossInfo.f(rep,:,:,:) = outRep.crossInfo.f;
%                 out.crossInfo.f1(rep,:,:,:) = outRep.crossInfo.f1;
%                 out.crossInfo.f2(rep,:,:,:) = outRep.crossInfo.f2;
                out.crossInfo_dim1Train = p.trainInfo;
%                 out.crossInfo_dim2Test = p.trainInfo;
%                 switch p.infoComb
%                     case 'layeredSVM'
%                         out.newRes.hitRate(rep,:) = outRep.newRes.hitRate;
%                         out.newRes.d(rep,:,:) = outRep.newRes.d;
%                         out.newRes.svmBias(rep,:,:) = outRep.newRes.svmBias;
%                         out.newRes.w(rep,:,:) = squeeze(outRep.newRes.w);
%                         out.newRes.a(rep,:,:) = squeeze(outRep.newRes.a);
%                         out.newRes.info(rep,:) = outRep.newRes.info;
%                     case 'layeredSVMwCross'
%                         out.combInfo.hitRate(rep,:) = outRep.combInfo.hitRate;
%                         out.combInfo.d(rep,:,:) = outRep.combInfo.decision';
%                     otherwise
%                         error('missing something');
%                 end
                clear outRep
                
%             case 'runSVM_patAct7'
% %                 keyboard;
% %                 p.trainInfo = {'p' 'm'}; p.testInfo = p.trainInfo;
%                 outRep = runSVM_patAct7(d,p,0,rep,doPerm);
%                 %Compile
%                 out.crossInfo.hitRate(rep,:,:,:) = outRep.crossInfo.hitRate;
%                 out.crossInfo.d(rep,:,:,:) = outRep.crossInfo.d;
%                 out.crossInfo.f(rep,:,:,:) = outRep.crossInfo.f;
% %                 out.crossInfo.f1(rep,:,:,:) = outRep.crossInfo.f1;
% %                 out.crossInfo.f2(rep,:,:,:) = outRep.crossInfo.f2;
%                 out.crossInfo_dim1Train = outRep.crossInfo_dim1Train;
%                 out.crossInfo_dim2Test = outRep.crossInfo_dim2Test;
%                 clear outRep
%             case 'runSVM_patAct8'
%                 error('not using this anymore, see within runSVM_patAct6.m')
% %                 
%                 %                 p.trainInfo = {'p'}; p.testInfo = {'r'};
%                 outRep = runSVM_patAct8(d,p,0,rep);
%                 %Compile
%                 out.crossInfo.hitRate(rep,:,:,:) = outRep.crossInfo.hitRate;
%                 out.crossInfo.d(rep,:,:,:) = outRep.crossInfo.d;
%                 out.crossInfo.w(rep,:,:,:) = outRep.crossInfo.w;
%                 out.crossInfo.a(rep,:,:,:) = outRep.crossInfo.a;
%                 out.crossInfo.f(rep,:,:,:) = outRep.crossInfo.f;
%                 out.crossInfo.f1(rep,:,:,:) = outRep.crossInfo.f1;
%                 out.crossInfo.f2(rep,:,:,:) = outRep.crossInfo.f2;
%                 out.crossInfo_dim1Train = outRep.crossInfo_dim1Train;
%                 out.crossInfo_dim2Test = outRep.crossInfo_dim2Test;
%                 clear outRep
        end
        
        if verbose
            dur = toc(time1);
            display(['Repetition ' num2str(rep) '; Took ' num2str(round(dur)) 'sec'])
        end
    end
end

switch p.algorithm
%     case 'runSVM_patAct'
%         svm.r.hitRate = out.hitRate;
%         svm.p = p;
%         svm.r.hitRate_d = out.hitRate_d;
%         svm.r.hitRate_pat = out.hitRate_pat;
%         svm.r.hitRate_patRel = out.hitRate_patRel;
%         svm.r.pat = out.pat;
%         svm.r.patRel = out.patRel;
    case {'runSVM_patAct2','runSVM_patAct4','runSVM_patAct5','runSVM_patAct6','runSVM_patAct7'}
        svm.r = out;
        svm.r.hitRate = [];
        svm.p = p;
    otherwise
        error('check that')
        svm.r.hitRate = hitRate;
        svm.p = p;
end


