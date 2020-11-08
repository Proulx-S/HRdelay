function svm = runSVMRepetitions3(d,p,verbose,usePara)
if ~exist('verbose','var')
    verbose = 1;
end
if ~exist('usePara','var')
    usePara = 0;
end


% if ischar(p.k)
%     switch p.k
%         case 'auto'
%             k = p.nObs/2;
%         case 'autoRun'
%             if strcmp(p.split(2:end),'perRun')
%                 k = p.nObs/2/str2double(p.split(1)); % to keep k constant (at the equivalent of leave-one-run-pair-out) across different spliting of the data
%             else
%                 error('')
%             end
%         otherwise
%             error('')
%     end
% else
%     k = p.k;
% end

% w = nan(p.nFeatLev,p.nFeatures,k,p.repeat);
% a = nan(p.nFeatLev,p.nFeatures,k,p.repeat);
hitRate = nan(p.repeat,p.nFeatLev,p.nFeatLev);
% if p.thresholdData; mi = nan(k,p.nFeatures,p.repeat); end
% amp = nan(p.nFeatLev,p.nFeatures,k,p.repeat);
% delay = nan(p.nFeatLev,p.nFeatures,k,p.repeat);

if usePara
    numworkers = 3;
    if verbose
        display(['Running SVM at ' num2str(length(p.featLevList)) '*' num2str(length(p.featLevList)) '/2 feature selection levels'])
        display(['for ' num2str(p.repeat) ' repetitions on ' num2str(numworkers) ' workers'])
    end
    
    poolobj = gcp('nocreate'); % If no pool, do not create new one.
    if isempty(poolobj)
        parpool(3);
    elseif poolobj.NumWorkers~=3
        delete(poolobj);
        parpool(3);
    else
        display(['parpool already created with ' num2str(numworkers) ' workers'])
    end
    
    parfor rep = 1:p.repeat
        %Run SVM
        hitRate(rep,:,:) = runSVM5(d,p,verbose,rep);
    end
    
else
    for rep = 1:p.repeat
        if verbose
            tic
            display(['Starting repetition ' num2str(rep) '/' num2str(p.repeat)])
        end
        
        %Run SVM
        hitRate(rep,:,:) = runSVM5(d,p,verbose,rep);
        
        if verbose
            dur = toc;
            display(['Repetition ' num2str(rep) '; Took ' num2str(round(dur)) 'sec'])
        end
    end
end



svm.r.hitRate = hitRate;
svm.p = p;

