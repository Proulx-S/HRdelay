function runPermDecoding(p)

permFlag = 1;
disp('***************************')
disp('Permutation test: computing')
disp('***************************')
for permInd = 1:p.perm.n
    disp(['Perm ' num2str(permInd) '/' num2str(p.perm.n)])
    if permInd==1
        resp = [];
        fov = [];
    end
    [respP,resp] = permuteLabels(p,resp);
    [featSelP,fov] = processFeatSel(p,permFlag,respP,fov);
    [resBSP,~,~,~,~] = runAllDecoding(p,permFlag,respP,featSelP);
    if permInd == 1
        auc = nan([size(resBSP) p.perm.n]);
    end
    tmp = reshape([resBSP{:}],size(resBSP));
%     tmp = reshape([tmp.subj],size(tmp));
    tmp = reshape(permute([tmp.auc],[2 3 1]),[size(tmp) size(tmp(1).auc,1)]);
    auc(:,:,permInd) = mean(tmp,3); clear tmp
end
disp('**********************')
disp('Permutation test: done')
disp('**********************')
% Save permutations
disp('Saving permutations')
auc = permute(auc,[3 1 2]);
resBSP = cell(size(resBSP));
for i = 1:numel(auc(1,:))
    resBSP{i}.subj.aucP = permute(auc(:,i),[2 3 1]);
end
outFile = fullfile(p.wd,'results',p.anaID,'decoding.mat');
save(outFile,'resBSP','-append');
disp(['Permutations appended to ' outFile ' as variable respBSP'])
