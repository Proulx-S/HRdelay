function svm = runSVMRepetitions5(d,p,verbose,ppn)
global tmpStrResp

if ~exist('verbose','var')
    verbose = 1;
end
if ~exist('ppn','var')
    ppn = 0;
end
% if ~exist('doPerm','var')
%     doPerm = 0;
%     p.doPerm = 0;
% else
%     p.doPerm = doPerm;
% end


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
        case 'paraFoldNsubFold'
            t1 = tic;
            display(['computing ' num2str(p.repeat) ' repetitions'])
            parfor rep = 1:p.repeat
                [rTr{rep},rVal{rep},rTe{rep}] = runSVMOneRepetition3(d,p,0,rep,tmpStrResp);
            end
            display(['computed in ' datestr(1/24/60/60*toc(t1), 'HH:MM:SS') 'sec'])
        otherwise
            error('X')
    end
    if isempty(getenv('OS'))
        delete(poolobj)
        rmdir(tmpJobStorageLocation,'s');
    end
else
    %% Run SVM with regular matlab
    % Run SVM
    switch p.algorithm
        case 'paraFoldNsubFold'
            for rep = 1:p.repeat
                if verbose
                    time1 = tic;
                    display(['Starting repetition ' num2str(rep) '/' num2str(p.repeat)])
                end
                
                [rTr{rep},rVal{rep},rTe{rep}] = runSVMOneRepetition3(d,p,0,rep);
                if isempty(rTr{rep})
                    svm = [];
                    return
                end
                
                if verbose
                    dur = toc(time1);
                    display(['Repetition ' num2str(rep) '; Took ' num2str(round(dur)) 'sec'])
                end
            end
        otherwise
            error('X')
    end
end

%% Compile repetitions
if isfield(rTr{1}{1},'cOpt')
    cOptFlag = 1;
    cOpt = rTr{1}{1}.cOpt;
    for rep = 1:length(rTr)
        for anaInd = 1:length(rTr{rep})
            rTr{rep}{anaInd} = rmfield(rTr{rep}{anaInd},'cOpt');
        end
    end
else
    cOptFlag = 0;
end

for anaInd = 1:length(rTr{1})
    if length(rTr{1})>1 && isfield(rTr{1}{anaInd},'pca')
        pcaFlag = 1;
        pcaRes{anaInd} = rTr{1}{anaInd}.pca;
        pcaRes{anaInd}.info = rTr{1}{anaInd}.info;
        for rep = 1:length(rTr)
            rTr{rep}{anaInd} = rmfield(rTr{rep}{anaInd},'pca');
        end
    else
        pcaFlag = 0;
    end
end
    
    
    
% if isfield(rTr{rep}{1},'pca')
%     pcaFlag = 1;
%     pcaRes = rTr{1}{1}.pca;
%     for rep = 1:length(rTr)
%         for anaInd = 1:length(rTr{rep})
%             rTr{rep}{anaInd} = rmfield(rTr{rep}{anaInd},'pca');
%         end
%     end
% else
%     pcaFlag = 0;
% end


if length(p.infoComb)==1 && strcmp(p.infoComb,'crossValDelay')
    for rep = 1:length(rTr)
        rTr2.diff.N(rep,:) = rTr{rep}{1}.diff.N;
        rTr2.diff.binCent(1,:) = rTr{rep}{1}.diff.binCent;
        rTr2.diff.binCent2(1,:) = rTr{rep}{1}.diff.binCent2;
        rTr2.diff.delay2(1,:) = rTr{rep}{1}.diff.delay2;
%         rTr2.t.N(rep,:) = rTr{rep}{1}.t.N;
%         rTr2.t.binCent(1,:) = rTr{rep}{1}.t.binCent;
        rTe2.diff.crossValDelay(rep,:) = rTe{rep}{1}.diff.crossValDelay;
        rTe2.diff.binCent(1,:) = rTe{rep}{1}.diff.binCent;
        rTe2.diff.crossValDelay2(rep,:) = rTe{rep}{1}.diff.crossValDelay2;
        rTe2.diff.binCent2(1,:) = rTe{rep}{1}.diff.binCent2;
%         rTe2.t.crossValDelay(rep,:) = rTe{rep}{1}.t.crossValDelay;
%         rTe2.t.binCent(1,:) = rTe{rep}{1}.t.binCent;
    end
    rVal2 = rVal;
elseif length(p.infoComb)==1 && strcmp(p.infoComb,'crossValDelayQ')
    for rep = 1:length(rTr)
        rTr2.diff.binCent(rep,:) = rTr{rep}{1}.diff.binCent;
        rTr2.diff.edges(rep,:) = rTr{rep}{1}.diff.edges;
        rTr2.diff.N(rep,:) = rTr{rep}{1}.diff.N;
        rTr2.diff.delay(rep,:) = rTr{rep}{1}.diff.delay;
        
        rTe2.diff.binCent(rep,:) = rTe{rep}{1}.diff.binCent;
        rTe2.diff.edges(rep,:) = rTe{rep}{1}.diff.edges;
        rTe2.diff.delay(rep,:) = rTe{rep}{1}.diff.delay;
        
    end
    rVal2 = rVal;
else
    rTr2 = rTr{1};
    if length(rTr)>1
        for rep = 2:length(rTr)
            rTr2 = compileRep2(rTr2,rTr{rep});
        end
    end
    rVal2 = rVal{1};
    if length(rVal)>1
        for rep = 2:length(rVal)
            rVal2 = compileRep2(rVal2,rVal{rep});
        end
    end
    rTe2 = rTe{1};
    if length(rTe)>1
        for rep = 2:length(rTr)
            rTe2 = compileRep2(rTe2,rTe{rep});
        end
    end
end



svm.d = d;
svm.p = p;
svm.rTr = rTr2;
svm.rVal = rVal2;
svm.rTe = rTe2;
if length(p.infoComb)==1 && strcmp(p.infoComb,'crossValDelay')
elseif length(p.infoComb)==1 && strcmp(p.infoComb,'crossValDelayQ')
else
    svm.rTe{1}.x = []; svm.rTe{1}.labelPairs = []; svm.rTe{1}.getPattern = [];
end
if cOptFlag
    svm.rTr{1}.cOpt = cOpt;
end
if pcaFlag
    for anaInd = 1:length(svm.rTr)
        svm.rTr{anaInd}.pca = pcaRes{anaInd};
    end
end



