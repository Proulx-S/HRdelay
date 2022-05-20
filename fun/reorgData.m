function [dTmp, info] = reorgData(p,d)
sessList = fields(d{1});
dTmp = cell(size(d,2),length(sessList));
for subjInd = 1:length(d)
    for sessInd = 1:length(sessList)
        dTmp{subjInd,sessInd} = d{subjInd}.(sessList{sessInd});
        d{subjInd}.(sessList{sessInd}) = [];
    end
end
info.info = 'subj x sess';
info.sessList = sessList';
info.subjList = p.meta.subjList;
