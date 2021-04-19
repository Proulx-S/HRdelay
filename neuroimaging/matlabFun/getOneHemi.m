function d = getOneHemi(d,hemi)

ind = d.voxProp.(['hemifield' hemi]);

fieldList = fields(d);
for fieldInd = 1:length(fieldList)
    sz = size(d.(fieldList{fieldInd}));
    if any(sz==length(ind))
        sz(sz==length(ind)) = nnz(ind);
        tmp = nan(sz);
        sz(sz==nnz(ind)) = 1;
        tmp(:) = d.(fieldList{fieldInd})(repmat(ind,sz));
        d.(fieldList{fieldInd}) = tmp; clear tmp sz
    end
end

fieldList = fields(d.voxProp);
for fieldInd = 1:length(fieldList)
    sz = size(d.voxProp.(fieldList{fieldInd}));
    if any(sz==length(ind))
        sz(sz==length(ind)) = nnz(ind);
        tmp = nan(sz);
        sz(sz==nnz(ind)) = 1;
        tmp(:) = d.voxProp.(fieldList{fieldInd})(repmat(ind,sz));
        d.voxProp.(fieldList{fieldInd}) = tmp; clear tmp sz
    end
end

d.voxProp.hemifieldL = logical(d.voxProp.hemifieldL);
d.voxProp.hemifieldR = logical(d.voxProp.hemifieldR);
d.voxProp.hemifield = ind;