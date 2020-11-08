function svm = runSVMRepetitions4(d,p,verbose,ppn,doPerm)
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
                    error('X');
            end
        else
            error('X')
            k = p.k;
            p.k = 'randK';
        end
        
        out.hitRate = nan(p.repeat,k);
        out.hitRate_d = nan(p.repeat,k);
        out.hitRate_pat = nan(p.repeat,k);
        out.hitRate_patRel = nan(p.repeat,k);
        out.pat = nan(p.repeat,p.nObs,2);
        out.patRel = nan(p.repeat,p.nObs,2);
    case {'runSVMOneRepetition','runSVMOneRepetition2'}
    otherwise
        error('X')
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
        case 'runSVM'
            parfor rep = 1:p.repeat
                hitRate(rep,:,:) = runSVM6(d,p,verbose,rep,doPerm);
            end
        case 'runSVM_RFE'
            parfor rep = 1:p.repeat
                hitRate(rep,:) = runSVM6(d,p,verbose,rep,doPerm);
            end
        case 'runSVM_patAct2'
            error('double-check that')
            parfor rep = 1:p.repeat
                outRep{rep} = runSVM_patAct3(d,p,0,rep,doPerm);
            end
            for rep = 1:p.repeat
                %Compile
                out.crossInfo.hitRate(rep,:,:) = outRep{rep}.crossInfo.hitRate;
                out.crossInfo.d(rep,:,:,:) = outRep{rep}.crossInfo.d;
                out.crossInfo_dim1Train = outRep{rep}.crossInfo_dim1Train;
                out.crossInfo_dim2Test = outRep{rep}.crossInfo_dim2Test;
            end
            clear outRep
        case 'runSVM_patAct4'
            parfor rep = 1:p.repeat
                outRep{rep} = runSVM_patAct4(d,p,0,rep,doPerm);
            end
            for rep = 1:p.repeat
                %Compile
                out.crossInfo.hitRate(rep,:,:,:) = outRep{rep}.crossInfo.hitRate;
                out.crossInfo.d(rep,:,:,:) = outRep{rep}.crossInfo.d;
                out.crossInfo.f(rep,:,:,:) = outRep{rep}.crossInfo.f;
                out.crossInfo.f1(rep,:,:,:) = outRep{rep}.crossInfo.f1;
                out.crossInfo.f2(rep,:,:,:) = outRep{rep}.crossInfo.f2;
                out.crossInfo_dim1Train = outRep{rep}.crossInfo_dim1Train;
                out.crossInfo_dim2Test = outRep{rep}.crossInfo_dim2Test;
            end
            clear outRep
        case 'runSVM_patAct5'
            parfor rep = 1:p.repeat
                outRep{rep} = runSVM_patAct5(d,p,0,rep,doPerm);
            end
            for rep = 1:p.repeat
                %Compile
                out.crossInfo.hitRate(rep,:,:,:) = outRep{rep}.crossInfo.hitRate;
                out.crossInfo.d(rep,:,:,:) = outRep{rep}.crossInfo.d;
                out.crossInfo.f(rep,:,:,:) = outRep{rep}.crossInfo.f;
                out.crossInfo.f1(rep,:,:,:) = outRep{rep}.crossInfo.f1;
                out.crossInfo.f2(rep,:,:,:) = outRep{rep}.crossInfo.f2;
                out.crossInfo_dim1Train = outRep{rep}.crossInfo_dim1Train;
                out.crossInfo_dim2Test = outRep{rep}.crossInfo_dim2Test;
            end
            clear outRep
        case 'runSVM_patAct6'
            parfor rep = 1:p.repeat
                outRep{rep} = runSVM_patAct6(d,p,0,rep);
            end
            for rep = 1:p.repeat
                %Compile
                out.crossInfo.hitRate(rep,:,:,:) = outRep{rep}.crossInfo.hitRate;
                out.crossInfo.d(rep,:,:,:) = outRep{rep}.crossInfo.d;
                out.crossInfo.w(rep,:,:,:) = outRep{rep}.crossInfo.w;
                out.crossInfo.a(rep,:,:,:) = outRep{rep}.crossInfo.a;
                out.crossInfo.f(rep,:,:,:) = outRep{rep}.crossInfo.f;
                out.crossInfo.f1(rep,:,:,:) = outRep{rep}.crossInfo.f1;
                out.crossInfo.f2(rep,:,:,:) = outRep{rep}.crossInfo.f2;
                out.crossInfo_dim1Train = outRep{rep}.crossInfo_dim1Train;
                out.crossInfo_dim2Test = outRep{rep}.crossInfo_dim2Test;
                switch p.infoComb
                    case 'layeredSVM'
                        out.newRes.hitRate(rep,:) = outRep{rep}.newRes.hitRate;
                        out.newRes.d(rep,:,:) = outRep{rep}.newRes.d;
                        out.newRes.svmBias(rep,:,:) = outRep{rep}.newRes.svmBias;
                        out.newRes.w(rep,:,:) = squeeze(outRep{rep}.newRes.w);
                        out.newRes.a(rep,:,:) = squeeze(outRep{rep}.newRes.a);
                        out.newRes.info(rep,:) = outRep{rep}.newRes.info;
                    case 'layeredSVMwCross'
                        out.combInfo.hitRate(rep,:) = outRep{rep}.combInfo.hitRate;
                        out.combInfo.d(rep,:,:) = outRep{rep}.combInfo.decision';
                    otherwise
                        error('missing something');
                end
            end
            clear outRep
        case 'runSVMOneRepetition'
            parfor rep = 1:p.repeat
                outRep{rep} = runSVMOneRepetition(d,p,0,rep);
            end
            for rep = 1:p.repeat
                %Compile
                out.crossInfo.hitRate(rep,:,:,:) = outRep{rep}.crossInfo.hitRate;
                out.crossInfo.d(rep,:,:,:) = outRep{rep}.crossInfo.d;
                out.crossInfo.w(rep,:,:,:) = permute(outRep{rep}.crossInfo.w,[1 3 2]);
                out.crossInfo.a(rep,:,:,:) = permute(outRep{rep}.crossInfo.a,[1 3 2]);
                out.crossInfo.f(rep,:,:,:) = outRep{rep}.crossInfo.f;
                out.crossInfo.f1(rep,:,:,:) = outRep{rep}.crossInfo.f1;
                out.crossInfo.f2(rep,:,:,:) = outRep{rep}.crossInfo.f2;
                out.crossInfo_dim1Train = outRep{rep}.crossInfo_dim1Train;
                out.crossInfo_dim2Test = outRep{rep}.crossInfo_dim2Test;
                if isfield(outRep{rep},'newRes')
                    out.newRes.hitRate(rep,:) = outRep{rep}.newRes.hitRate;
                    out.newRes.d(rep,:,:) = outRep{rep}.newRes.d;
                    for i = 1:size(outRep{rep}.newRes.w,1)
                        out.newRes.w{i}(rep,:) = outRep{rep}.newRes.w{i};
                        out.newRes.a{i}(rep,:) = outRep{rep}.newRes.a{i};
                    end
                    
                    out.newRes.hitRateL2O1(rep,:) = outRep{rep}.newRes.hitRateL2O1;
                    out.newRes.dL2O1(rep,:,:) = outRep{rep}.newRes.dL2O1;
                    for i = 1:size(outRep{rep}.newRes.wL2O1,1)
                        if isempty(outRep{rep}.newRes.wL2O1{i})
                            out.newRes.wL2O1{i} = 'n/a';
                            out.newRes.aL2O1{i} = 'n/a';
                        else
                            out.newRes.wL2O1{i}(rep,:) = outRep{rep}.newRes.wL2O1{i};
                            out.newRes.aL2O1{i}(rep,:) = outRep{rep}.newRes.aL2O1{i};
                        end
                    end
                    
                    out.newRes.hitRateL2O2(rep,:) = outRep{rep}.newRes.hitRateL2O2;
                    out.newRes.dL2O2(rep,:,:) = outRep{rep}.newRes.dL2O2;
                    for i = 1:size(outRep{rep}.newRes.wL2O2,1)
                        if isempty(outRep{rep}.newRes.wL2O2{i})
                            out.newRes.wL2O2{i} = 'n/a';
                            out.newRes.aL2O2{i} = 'n/a';
                        else
                            out.newRes.wL2O2{i}(rep,:) = outRep{rep}.newRes.wL2O2{i};
                            out.newRes.aL2O2{i}(rep,:) = outRep{rep}.newRes.aL2O2{i};
                        end
                    end
                    
                    out.newRes.hitRateL1O1(rep,:) = outRep{rep}.newRes.hitRateL1O1;
                    out.newRes.dL1O1(rep,:,:) = outRep{rep}.newRes.dL1O1;
                    for i = 1:size(outRep{rep}.newRes.wL1O1,1)
                        if isempty(outRep{rep}.newRes.wL1O1{i})
                            out.newRes.wL1O1{i} = 'n/a';
                            out.newRes.aL1O1{i} = 'n/a';
                        else
                            out.newRes.wL1O1{i}(rep,:) = outRep{rep}.newRes.wL1O1{i};
                            out.newRes.aL1O1{i}(rep,:) = outRep{rep}.newRes.aL1O1{i};
                        end
                    end
                    
                    out.newRes.hitRateL1O2(rep,:) = outRep{rep}.newRes.hitRateL1O2;
                    out.newRes.dL1O2(rep,:,:) = outRep{rep}.newRes.dL1O2;
                    for i = 1:size(outRep{rep}.newRes.wL1O2,1)
                        if isempty(outRep{rep}.newRes.wL1O2{i})
                            out.newRes.wL1O2{i} = 'n/a';
                            out.newRes.aL1O2{i} = 'n/a';
                        else
                            out.newRes.wL1O2{i}(rep,:) = outRep{rep}.newRes.wL1O2{i};
                            out.newRes.aL1O2{i}(rep,:) = outRep{rep}.newRes.aL1O2{i};
                        end
                    end
                    
                    %                     out.newRes.svmBias(rep,:,:) = outRep{rep}.newRes.svmBias;
                    out.newRes.info(rep,:) = outRep{rep}.newRes.info;
                    out.newRes.infoL2O1(rep,:) = outRep{rep}.newRes.infoL2O1;
                    out.newRes.infoL2O2(rep,:) = outRep{rep}.newRes.infoL2O2;
                    out.newRes.infoL1O1(rep,:) = outRep{rep}.newRes.infoL1O1;
                    out.newRes.infoL1O2(rep,:) = outRep{rep}.newRes.infoL1O2;
                end
            end
            clear outRep
        otherwise
            error('X')
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
            case 'runSVM'
                hitRate(rep,:,:) = runSVM6(d,p,verbose,rep,doPerm);
            case 'runSVM_RFE'
                keyboard
                hitRate(rep,:) = runSVM6(d,p,verbose,rep,doPerm);
            case 'runSVM_patAct'
                outRep = runSVM_patAct2(d,p,verbose,rep,doPerm);
                %Compile
                allFields = fields(outRep);
                for i = 1:length(allFields)
                    if strfind(allFields{i},'hitRate')
                        out.(allFields{i})(rep,:) = outRep.(allFields{i});
                    elseif strfind(allFields{i},'pat')
                        out.(allFields{i})(rep,:,:) = outRep.(allFields{i});
                    else
                        error('something wrong here')
                    end
                end
                clear outRep
            case 'runSVM_patAct2'
                error('double-check that')
                outRep = runSVM_patAct3(d,p,0,rep,doPerm);
                
                %Compile
                out.crossInfo.hitRate(rep,:,:) = outRep.crossInfo.hitRate;
                out.crossInfo.d(rep,:,:,:) = outRep.crossInfo.d;
                out.crossInfo_dim1Train = outRep.crossInfo_dim1Train;
                out.crossInfo_dim2Test = outRep.crossInfo_dim2Test;
                clear outRep
            case 'runSVM_patAct4'
                outRep = runSVM_patAct4(d,p,0,rep,doPerm);
                %Compile
                out.crossInfo.hitRate(rep,:,:,:) = outRep.crossInfo.hitRate;
                out.crossInfo.d(rep,:,:,:) = outRep.crossInfo.d;
                out.crossInfo.f(rep,:,:,:) = outRep.crossInfo.f;
                out.crossInfo.f1(rep,:,:,:) = outRep.crossInfo.f1;
                out.crossInfo.f2(rep,:,:,:) = outRep.crossInfo.f2;
                out.crossInfo_dim1Train = outRep.crossInfo_dim1Train;
                out.crossInfo_dim2Test = outRep.crossInfo_dim2Test;
                clear outRep
            case 'runSVM_patAct5'
                outRep = runSVM_patAct5(d,p,0,rep,doPerm);
                %Compile
                out.crossInfo.hitRate(rep,:,:,:) = outRep.crossInfo.hitRate;
                out.crossInfo.d(rep,:,:,:) = outRep.crossInfo.d;
                out.crossInfo.f(rep,:,:,:) = outRep.crossInfo.f;
                out.crossInfo.f1(rep,:,:,:) = outRep.crossInfo.f1;
                out.crossInfo.f2(rep,:,:,:) = outRep.crossInfo.f2;
                out.crossInfo_dim1Train = outRep.crossInfo_dim1Train;
                out.crossInfo_dim2Test = outRep.crossInfo_dim2Test;
                clear outRep
            case 'runSVM_patAct6'
                %                 p.trainInfo = {'p'}; p.testInfo = {'r'};
                outRep = runSVM_patAct6(d,p,0,rep);
                %Compile
                out.crossInfo.hitRate(rep,:,:,:) = outRep.crossInfo.hitRate;
                out.crossInfo.d(rep,:,:,:) = outRep.crossInfo.d;
                out.crossInfo.w(rep,:,:,:) = outRep.crossInfo.w;
                out.crossInfo.a(rep,:,:,:) = outRep.crossInfo.a;
                out.crossInfo.f(rep,:,:,:) = outRep.crossInfo.f;
                out.crossInfo.f1(rep,:,:,:) = outRep.crossInfo.f1;
                out.crossInfo.f2(rep,:,:,:) = outRep.crossInfo.f2;
                out.crossInfo_dim1Train = outRep.crossInfo_dim1Train;
                out.crossInfo_dim2Test = outRep.crossInfo_dim2Test;
                switch p.infoComb
                    case 'layeredSVM'
                        out.newRes.hitRate(rep,:) = outRep.newRes.hitRate;
                        out.newRes.d(rep,:,:) = outRep.newRes.d;
                        out.newRes.svmBias(rep,:,:) = outRep.newRes.svmBias;
                        out.newRes.w(rep,:,:) = squeeze(outRep.newRes.w);
                        out.newRes.a(rep,:,:) = squeeze(outRep.newRes.a);
                        out.newRes.info(rep,:) = outRep.newRes.info;
                    case 'layeredSVMwCross'
                        out.combInfo.hitRate(rep,:) = outRep.combInfo.hitRate;
                        out.combInfo.d(rep,:,:) = outRep.combInfo.decision';
                    otherwise
                        error('missing something');
                end
                clear outRep
                
            case 'runSVMOneRepetition2'
                outRep = runSVMOneRepetition2(d,p,0,rep);
            case 'runSVMOneRepetition'
                %                 p.trainInfo = {'p'}; p.testInfo = {'r'};
                outRep = runSVMOneRepetition(d,p,0,rep);
                %Compile
                out.crossInfo.hitRate(rep,:,:,:) = outRep.crossInfo.hitRate;
                out.crossInfo.d(rep,:,:,:) = outRep.crossInfo.d;
                out.crossInfo.w(rep,:,:,:) = permute(outRep.crossInfo.w,[1 3 2]);
                out.crossInfo.a(rep,:,:,:) = permute(outRep.crossInfo.a,[1 3 2]);
                out.crossInfo.f(rep,:,:,:) = outRep.crossInfo.f;
                out.crossInfo.f1(rep,:,:,:) = outRep.crossInfo.f1;
                out.crossInfo.f2(rep,:,:,:) = outRep.crossInfo.f2;
                out.crossInfo_dim1Train = outRep.crossInfo_dim1Train;
                out.crossInfo_dim2Test = outRep.crossInfo_dim2Test;
                if isfield(outRep,'newRes')
                    out.newRes.hitRate(rep,:) = outRep.newRes.hitRate;
                    out.newRes.d(rep,:,:) = outRep.newRes.d;
                    for i = 1:size(outRep.newRes.w,1)
                        out.newRes.w{i}(rep,:) = outRep.newRes.w{i};
                        out.newRes.a{i}(rep,:) = outRep.newRes.a{i};
                    end
                    
                    out.newRes.hitRateL2O1(rep,:) = outRep.newRes.hitRateL2O1;
                    out.newRes.dL2O1(rep,:,:) = outRep.newRes.dL2O1;
                    for i = 1:size(outRep.newRes.wL2O1,1)
                        if isempty(outRep.newRes.wL2O1{i})
                            out.newRes.wL2O1{i} = 'n/a';
                            out.newRes.aL2O1{i} = 'n/a';
                        else
                            out.newRes.wL2O1{i}(rep,:) = outRep.newRes.wL2O1{i};
                            out.newRes.aL2O1{i}(rep,:) = outRep.newRes.aL2O1{i};
                        end
                    end
                    
                    out.newRes.hitRateL2O2(rep,:) = outRep.newRes.hitRateL2O2;
                    out.newRes.dL2O2(rep,:,:) = outRep.newRes.dL2O2;
                    for i = 1:size(outRep.newRes.wL2O2,1)
                        if isempty(outRep.newRes.wL2O2{i})
                            out.newRes.wL2O2{i} = 'n/a';
                            out.newRes.aL2O2{i} = 'n/a';
                        else
                            out.newRes.wL2O2{i}(rep,:) = outRep.newRes.wL2O2{i};
                            out.newRes.aL2O2{i}(rep,:) = outRep.newRes.aL2O2{i};
                        end
                    end
                    
                    out.newRes.hitRateL1O1(rep,:) = outRep.newRes.hitRateL1O1;
                    out.newRes.dL1O1(rep,:,:) = outRep.newRes.dL1O1;
                    for i = 1:size(outRep.newRes.wL1O1,1)
                        if isempty(outRep.newRes.wL1O1{i})
                            out.newRes.wL1O1{i} = 'n/a';
                            out.newRes.aL1O1{i} = 'n/a';
                        else
                            out.newRes.wL1O1{i}(rep,:) = outRep.newRes.wL1O1{i};
                            out.newRes.aL1O1{i}(rep,:) = outRep.newRes.aL1O1{i};
                        end
                    end
                    
                    out.newRes.hitRateL1O2(rep,:) = outRep.newRes.hitRateL1O2;
                    out.newRes.dL1O2(rep,:,:) = outRep.newRes.dL1O2;
                    for i = 1:size(outRep.newRes.wL1O2,1)
                        if isempty(outRep.newRes.wL1O2{i})
                            out.newRes.wL1O2{i} = 'n/a';
                            out.newRes.aL1O2{i} = 'n/a';
                        else
                            out.newRes.wL1O2{i}(rep,:) = outRep.newRes.wL1O2{i};
                            out.newRes.aL1O2{i}(rep,:) = outRep.newRes.aL1O2{i};
                        end
                    end
                    
                    %                     out.newRes.svmBias(rep,:,:) = outRep.newRes.svmBias;
                    out.newRes.info(rep,:) = outRep.newRes.info;
                    out.newRes.infoL2O1(rep,:) = outRep.newRes.infoL2O1;
                    out.newRes.infoL2O2(rep,:) = outRep.newRes.infoL2O2;
                    out.newRes.infoL1O1(rep,:) = outRep.newRes.infoL1O1;
                    out.newRes.infoL1O2(rep,:) = outRep.newRes.infoL1O2;
                end
                clear outRep
                
            case 'runSVM_patAct7'
                %                 keyboard;
                %                 p.trainInfo = {'p' 'm'}; p.testInfo = p.trainInfo;
                outRep = runSVM_patAct7(d,p,0,rep,doPerm);
                %Compile
                out.crossInfo.hitRate(rep,:,:,:) = outRep.crossInfo.hitRate;
                out.crossInfo.d(rep,:,:,:) = outRep.crossInfo.d;
                out.crossInfo.f(rep,:,:,:) = outRep.crossInfo.f;
                %                 out.crossInfo.f1(rep,:,:,:) = outRep.crossInfo.f1;
                %                 out.crossInfo.f2(rep,:,:,:) = outRep.crossInfo.f2;
                out.crossInfo_dim1Train = outRep.crossInfo_dim1Train;
                out.crossInfo_dim2Test = outRep.crossInfo_dim2Test;
                clear outRep
            case 'runSVM_patAct8'
                error('not using this anymore, see within runSVM_patAct6.m')
                %
                %                 p.trainInfo = {'p'}; p.testInfo = {'r'};
                outRep = runSVM_patAct8(d,p,0,rep);
                %Compile
                out.crossInfo.hitRate(rep,:,:,:) = outRep.crossInfo.hitRate;
                out.crossInfo.d(rep,:,:,:) = outRep.crossInfo.d;
                out.crossInfo.w(rep,:,:,:) = outRep.crossInfo.w;
                out.crossInfo.a(rep,:,:,:) = outRep.crossInfo.a;
                out.crossInfo.f(rep,:,:,:) = outRep.crossInfo.f;
                out.crossInfo.f1(rep,:,:,:) = outRep.crossInfo.f1;
                out.crossInfo.f2(rep,:,:,:) = outRep.crossInfo.f2;
                out.crossInfo_dim1Train = outRep.crossInfo_dim1Train;
                out.crossInfo_dim2Test = outRep.crossInfo_dim2Test;
                clear outRep
            otherwise
                error('X')
        end
        
        if verbose
            dur = toc(time1);
            display(['Repetition ' num2str(rep) '; Took ' num2str(round(dur)) 'sec'])
        end
    end
    xxx=1;
end

switch p.algorithm
    case 'runSVM_patAct'
        svm.r.hitRate = out.hitRate;
        svm.p = p;
        svm.r.hitRate_d = out.hitRate_d;
        svm.r.hitRate_pat = out.hitRate_pat;
        svm.r.hitRate_patRel = out.hitRate_patRel;
        svm.r.pat = out.pat;
        svm.r.patRel = out.patRel;
    case {'runSVM_patAct2','runSVM_patAct4','runSVM_patAct5','runSVM_patAct6','runSVM_patAct7','runSVMOneRepetition'}
        svm.r = out;
        svm.r.hitRate = out.crossInfo.hitRate;
        svm.p = p;
    otherwise
        error('check that')
        svm.r.hitRate = hitRate;
        svm.p = p;
end


