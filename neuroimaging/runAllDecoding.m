function [resBS,resBShr,resWS,f,info] = runAllDecoding(p,verbose)
if ~exist('verbose','var')
    verbose = 1;
end
if ~isfield(p,'perm')
    p.perm.doIt = 0;
end


%% Define paths
subjList = p.meta.subjList;
repoPath = p.paths.repo.in;
    funPath = fullfile(repoPath,'C-derived\DecodingHR\fun');
        outDir  = ['e_' p.anaID];
%make sure everything is forward slash for mac, linux pc compatibility
for tmp = {'repoPath' 'funPath' 'outDir'}
    eval([char(tmp) '(strfind(' char(tmp) ',''\''))=''/'';']);
end
clear tmp




%%
figOption = p.figOption;
condPairList = p.svm.condPairList;
respFeatList = p.svm.respFeatList;



figOption_verbose = p.figOption.verbose;


resBS = cell(length(condPairList),length(respFeatList));
resBShr = cell(length(condPairList),length(respFeatList));
resWS = cell(length(condPairList),length(respFeatList));
f = cell(length(condPairList),length(respFeatList));
p.complexSpace = p.svm.complexSpace;
for respFeatInd = 1:length(respFeatList)
    p.chanSpace = respFeatList{respFeatInd};
    for condPairInd = 1:length(condPairList)
        if condPairInd~=1
            p.figOption.verbose = 0;
        else
            p.figOption.verbose = figOption_verbose;
        end
        p.condPair = condPairList{condPairInd};
        [resBS{condPairInd,respFeatInd},resBShr{condPairInd,respFeatInd},resWS{condPairInd,respFeatInd},f{condPairInd,respFeatInd}] = runDecoding(p,verbose);
    end
end
p.figOption.verbose = figOption_verbose;
info.condPairList = condPairList';
info.respFeatList = respFeatList;
info.info = 'condPair x respFeat';

if ~p.perm.doIt
    %% Save data
    fullpath = fullfile(funPath,outDir);
    if ~exist(fullpath,'dir'); mkdir(fullpath); end
    fullfilename = fullfile(fullpath,'decoding');
    save(fullfilename,'resBS','f','p','info');
    fullfilename = fullfile(fullpath,'channels');
    save(fullfilename,'resBShr','f','p','info');
    
    %% Save figures
    if p.figOption.verbose
        f = mat2cell([f{:}],1,ones(1,numel([f{:}])));
        fullfilename = fullfile(p.figOption.finalDir,'decodingFig');
        for i = 1:numel(f)
            curF = f{i};
            curF.Color = 'none';
            set(findobj(curF.Children,'type','Axes'),'color','none')
            curFile = [fullfilename '_' num2str(i)];
            curExt = 'svg';
            saveas(curF,[curFile '.' curExt]); if p.figOption.verbose; disp([curFile '.' curExt]); end
            curF.Color = 'w';
            curExt = 'fig';
            saveas(curF,[curFile '.' curExt]); if p.figOption.verbose; disp([curFile '.' curExt]); end
            curExt = 'jpg';
            saveas(curF,[curFile '.' curExt]); if p.figOption.verbose; disp([curFile '.' curExt]); end
        end
    end
end

