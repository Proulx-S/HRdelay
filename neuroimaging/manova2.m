function [anovatbl_MV,A,C_MV,D] = manova2(rm_WM,withinModel,By,)

%% Run standard matlab manova
if exist('By','var') && ~isempty(By)
    ByFlag = 1;
    [anovatbl,A,C,D] = manova(rm_WM,'WithinModel',withinModel,'By',By);
else
    ByFlag = 0;
    [anovatbl,A,C,D] = manova(rm_WM,'WithinModel',withinModel);
end
% [anovatbl,A,C,D] = manova(rm_WM);
% Keep only Hotelling
anovatbl = anovatbl(3:4:end,:);
%% Redefine contrasts with the dependent variable factor (the first in withinModel_MV)
[C_MV,C_MV_names] = redefineC(rm_WM,anovatbl,C);

%% Recalculate stats
for betweenTestInd = 1:length(A)
    A_names{betweenTestInd} = char(table2array(anovatbl(betweenTestInd,2)));
end
withinNames = [];
betweenNames = [];
anovatbl_MV = [];
for withinTestInd = 1:length(C_MV)
    for betweenTestInd = 1:length(A)
        tmp = coeftest(rm_WM,A{betweenTestInd},C_MV{withinTestInd},D);
        anovatbl_MV = cat(1,anovatbl_MV,tmp(3,:));
%         withinNames = cat(1,withinNames,repmat(C_MV_names(withinTestInd),[4 1]));
        withinNames = cat(1,withinNames,C_MV_names(withinTestInd));
        betweenNames = cat(1,betweenNames,A_names(betweenTestInd));
%         if betweenTestInd==1
%             if ByFlag
%                 betweenNames = cat(1,betweenNames,repmat({char(table2array(anovatbl(1,2)))},[4 1]));
%             else
%                 betweenNames = cat(1,betweenNames,repmat({'(Intercept)'},[4 1]));
%             end
%         else
%             if ByFlag
%                 betweenNames = cat(1,betweenNames,repmat({char(table2array(anovatbl(5,2)))},[4 1]));
%             else
%                 betweenNames = cat(1,betweenNames,repmat(rm_WM.BetweenFactorNames(betweenTestInd-1),[4 1]));
%             end
%         end
    end 
end
tmp = table(withinNames,betweenNames,'VariableNames',anovatbl(:,1:2).Properties.VariableNames);
anovatbl_MV = cat(2,tmp,anovatbl_MV);




function [C,C_names] = redefineC(rm,ranovatbl,C)

depVarName = rm.WithinFactorNames{1};

switch length(rm.BetweenFactorNames)
    case 0
        btwTerms = 1;
    case 1
        btwTerms = length(rm.BetweenFactorNames)+1;
    case 2
        btwTerms = length(rm.BetweenFactorNames)+1+1;
    otherwise
        error('did not properly code for more than 2 between-subject factors')
end
testNames = table2cell(ranovatbl(1:btwTerms:end,1)); for i = 1:length(testNames); testNames(i) = cellstr(char(testNames{i})); end
CmulInd = find(~cellfun('isempty',strfind(testNames,depVarName)));
CmulNames = testNames(CmulInd);
for i = 1:length(CmulNames)
    curName = CmulNames{i};
    if ~contains(curName,[depVarName ':'])
        curNameInd = CmulInd(i);
        curNameToAdd = '(Intercept)';
        curNameToAddInd = find(ismember(testNames,curNameToAdd));
        if length(curNameToAddInd)>1
            error('X')
        end
        C{curNameInd} = [C{curNameToAddInd} C{curNameInd}];
        CmulNames{i} = '(Intercept)';
    else
        curNameInd = CmulInd(i);
        curNameToAdd = strrep(curName,[depVarName ':'],'');
        curNameToAddInd = find(ismember(testNames,curNameToAdd));
        if length(curNameToAddInd)>1
            error('X')
        end
        C{curNameInd} = [C{curNameToAddInd} C{curNameInd}];
        CmulNames{i} = curNameToAdd;
    end
end
C = C(CmulInd);
C_names = CmulNames';



