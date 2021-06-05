function visualizeOthers(p,verbose)
if ~exist('verbose','var')
    verbose = 1;
end
% if ~exist('figOption','var') || isempty(figOption)
%     figOption.save = 0;
%     figOption.subj = 1; % 'all' or subjInd
% end
if ~isfield(p,'chanSpace') || isempty(p.chanSpace)
    p.chanSpace = 'cart'; % 'cart_HTbSess' 'cartNoAmp_HTbSess' 'cartNoDelay_HTbSess'
    % 'hr' 'hrNoAmp' 'cart' 'cartNoAmp' cartNoAmp_HT 'cartReal', 'cartImag', 'pol', 'polMag' 'polMag_T' or 'polDelay'
end
f = [];
% tic

%% Define paths
subjList = {'02jp' '03sk' '04sp' '05bm' '06sb' '07bj'};
if ismac
    repoPath = '/Users/sebastienproulx/OneDrive - McGill University/dataBig';
else
    repoPath = 'C:\Users\sebas\OneDrive - McGill University\dataBig';
end
funPath = fullfile(repoPath,'C-derived\DecodingHR\fun');
inDir  = 'd';
%make sure everything is forward slash for mac, linux pc compatibility
for tmp = {'repoPath' 'funPath' 'inDir'}
    eval([char(tmp) '(strfind(' char(tmp) ',''\''))=''/'';']);
end
clear tmp



%% Load data
dAll = cell(size(subjList,1),1);
for subjInd = 1:size(subjList,2)
    curFile = fullfile(funPath,inDir,[subjList{subjInd} '.mat']);
    if verbose; disp(['loading: ' curFile]); end
    load(curFile,'res');
    dAll{subjInd} = res;
end
d = dAll; clear dAll
sessList = fields(d{1});
% Load feature slection
load(fullfile(funPath,inDir,'featSel.mat'),'featSel');
% if verbose
%     disp('---');
%     disp(['Channel space: ' p.chanSpace '-' p.condPair]);
%     disp(['Complex space: ' p.complexSpace]);
% end

%% Reorganize
dP = cell(size(d,2),length(sessList));
for subjInd = 1:length(d)
    for sessInd = 1:length(sessList)
        dP{subjInd,sessInd} = d{subjInd}.(sessList{sessInd});
        d{subjInd}.(sessList{sessInd}) = [];
        dP{subjInd,sessInd}.featSel = featSel{subjInd,sessInd};
    end
end
d = dP; clear dP

% line_num=dbstack; % get this line number
% disp(['line ' num2str(line_num(1).line)]) % displays the line number
% toc

%% Define feature selection
if ~strcmp(featSel{1,1}.featSeq.info2,p.featSel.global.method)
    error(['p.featSel.global.method and featSel.featSeq.info2 not matching' newline 'Try reruning processFeatSel.m'])
end
featSelSteps_labelList = featSel{1,1}.featSeq.featSelList;
featSelConds_labelList = featSel{1,1}.featSeq.condPairList;


method = 'onlyRetinoFov'; % 'onlyRetinoFov' 'upToActivation'
condPair = 'all';
[ind_nSpecFeatSel,ind_nSpecFeatSelCond,ind_specFeatSel,ind_specFeatSelCond] = defineFeatSel(featSelSteps_labelList,featSelConds_labelList,method,condPair);
% [ind_nSpecFeatSel,ind_nSpecFeatSelCond,ind_specFeatSel,ind_specFeatSelCond] = defineFeatSel(featSelSteps_labelList,featSelConds_labelList,p.featSel.global.method,p.condPair);
for subjInd = 1:size(d,1)
    for sessInd = 1:size(d,2)
        featSel{subjInd,sessInd}.indIn = ...
            all(featSel{subjInd,sessInd}.featSeq.featIndIn(:,ind_nSpecFeatSel,ind_nSpecFeatSelCond),2)...
            & all(featSel{subjInd,sessInd}.featSeq.featIndIn(:,ind_specFeatSel,ind_specFeatSelCond),2);
    end
end
featSel_fov = featSel;


method = 'upToActivation'; % 'onlyRetinoFov' 'upToActivation'
condPair = 'all';
[ind_nSpecFeatSel,ind_nSpecFeatSelCond,ind_specFeatSel,ind_specFeatSelCond] = defineFeatSel(featSelSteps_labelList,featSelConds_labelList,method,condPair);
% [ind_nSpecFeatSel,ind_nSpecFeatSelCond,ind_specFeatSel,ind_specFeatSelCond] = defineFeatSel(featSelSteps_labelList,featSelConds_labelList,p.featSel.global.method,p.condPair);
for subjInd = 1:size(d,1)
    for sessInd = 1:size(d,2)
        featSel{subjInd,sessInd}.indIn = ...
            all(featSel{subjInd,sessInd}.featSeq.featIndIn(:,ind_nSpecFeatSel,ind_nSpecFeatSelCond),2)...
            & all(featSel{subjInd,sessInd}.featSeq.featIndIn(:,ind_specFeatSel,ind_specFeatSelCond),2);
    end
end
featSel_act = featSel;

% line_num=dbstack; % get this line number
% disp(['line ' num2str(line_num(1).line)]) % displays the line number
% toc

%% Example plot of trigonometric (polar) representation
subjInd = p.figOption.subjInd;
sessInd = p.figOption.sessInd;
if p.figOption.verbose==1
    f = [f plotTrig_featSel(d{subjInd,sessInd},p,featSel_act{subjInd,sessInd},featSel_fov{subjInd,sessInd}),0];
elseif p.figOption.verbose>1
    f = [f plotTrig_featSel(d{subjInd,sessInd},p,featSel_act{subjInd,sessInd},featSel_fov{subjInd,sessInd}),1];
end
% if p.figOption.verbose>=1 && figOption.save
%     error('code that')
%     filename = fullfile(pwd,mfilename);
%     if ~exist(filename,'dir'); mkdir(filename); end
%     filename = fullfile(filename,p.chanSpace);
%     f.Color = 'none';
%     set(findobj(f.Children,'type','Axes'),'color','none')
%     set(findobj(f.Children,'type','PolarAxes'),'color','none')
%     saveas(f,[filename '.svg']); if verbose; disp([filename '.svg']); end
%     f.Color = 'w';
%     set(findobj(f.Children,'type','Axes'),'color','w')
%     set(findobj(f.Children,'type','PolarAxes'),'color','w')
%     saveas(f,[filename '.fig']); if verbose; disp([filename '.fig']); end
%     saveas(f,[filename '.jpg']); if verbose; disp([filename '.jpg']); end
% end

% line_num=dbstack; % get this line number
% disp(['line ' num2str(line_num(1).line)]) % displays the line number
% toc


f = [f plotNorm(d{subjInd,sessInd},p,featSel{subjInd,sessInd},[],0)];

