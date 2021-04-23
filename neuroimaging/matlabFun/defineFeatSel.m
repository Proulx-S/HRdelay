function [ind_nSpecFeatSel,ind_nSpecFeatSelCond,ind_specFeatSel,ind_specFeatSelCond] = defineFeatSel(featSel,p)


featSelSteps_labelList = featSel{1,1}.featSeq.featSelList;
for i = 1:length(featSelSteps_labelList)
    tmp = strsplit(featSelSteps_labelList{i},':');
    featSelSteps_labelList(i) = tmp(1);
end
featSelConds_labelList = featSel{1,1}.featSeq.condPairList;
for i = 1:length(featSelConds_labelList)
    featSelConds_labelList{i} = num2str(featSelConds_labelList{i});
end
switch p.featSel.global.method
    case 'allData'
        % Non-condition-specific featSel steps
        ind_nSpecFeatSel = ismember(featSelSteps_labelList,{'retinoFov' 'act' 'respVecSig' 'vein' 'respVecDiff'});
        % which condition set
        ind_nSpecFeatSelCond = [1 2 3];
        % Condition-specific fetSel steps
        ind_specFeatSel = ismember(featSelSteps_labelList,{''});
    case 'custom1'
        % Non-condition-specific featSel steps
        ind_nSpecFeatSel = ismember(featSelSteps_labelList,{'retinoFov' 'act' 'respVecSig' 'vein'});
        % which condition set
        ind_nSpecFeatSelCond = [1 2 3];
        % Condition-specific fetSel steps
        ind_specFeatSel = ismember(featSelSteps_labelList,{'respVecDiff'});
        % which condition set
        switch p.condPair
            case 'grat1VSgrat2'
                ind_specFeatSelCond = [1 2];
            case 'grat1VSplaid'
                ind_specFeatSelCond = [1 3];
            case 'grat2VSplaid'
                ind_specFeatSelCond = [2 3];
            otherwise
                error('X')
        end
    case 'custom2'
        % Non-condition-specific featSel steps
        ind_nSpecFeatSel = ismember(featSelSteps_labelList,{'retinoFov' 'vein'});
        % which condition set
        ind_nSpecFeatSelCond = [1 2 3];
        % Condition-specific fetSel steps
        ind_specFeatSel = ismember(featSelSteps_labelList,{'act' 'respVecSig' 'respVecDiff'});
        % which condition set
        switch p.condPair
            case 'grat1VSgrat2'
                ind_specFeatSelCond = [1 2];
            case 'grat1VSplaid'
                ind_specFeatSelCond = [1 3];
            case 'grat2VSplaid'
                ind_specFeatSelCond = [2 3];
            otherwise
                error('X')
        end
    otherwise
        error('X')
end
ind_nSpecFeatSelCond = squeeze(ismember(featSelConds_labelList,num2str(ind_nSpecFeatSelCond)));
ind_specFeatSelCond = squeeze(ismember(featSelConds_labelList,num2str(ind_specFeatSelCond)));
