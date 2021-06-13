function permuteLabels(p)
repoPath = p.paths.repo.in;
        funPath = fullfile(repoPath,'C-derived\DecodingHR\fun');
            inDir  = 'd';
            outDir = 'd';
subjList = p.meta.subjList;
sessList = {'sess1' 'sess2'};

disp('Permuting Labels: Permuting')
for subjInd = 1:length(subjList)
    subj = subjList{subjInd};
    %Load
    load(fullfile(funPath,inDir,[subj '.mat']),'res')
    for sessInd = 1:length(sessList)
        sess = sessList{sessInd};
        resP.(sess).sin = res.(sess).sin;
        resP.(sess).sinBase = res.(sess).sinBase;
        resP.(sess).hr = res.(sess).hr;
        resP.(sess).hrBase = res.(sess).hrBase;
        for rep = 1:size(res.(sess).sin,4)
            %Permute
            perm = randperm(3);
            resP.(sess).sin(:,:,:,rep,:) = res.(sess).sin(:,:,:,rep,perm);
            resP.(sess).sinBase(:,:,:,rep,:) = res.(sess).sinBase(:,:,:,rep,perm);
            resP.(sess).hr(:,:,:,rep,:,:) = res.(sess).hr(:,:,:,rep,perm,:);
            resP.(sess).hrBase(:,:,:,rep,:,:) = res.(sess).hrBase(:,:,:,rep,perm,:);
        end
    end
    %Save
    save(fullfile(funPath,outDir,[subj '.mat']),'resP','-append')
end
disp('Permuting Labels: Done')