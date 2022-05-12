function p = getRepoData(p)

if isempty(p.paths.data)
    [a,b,~] = fileparts(p.paths.pwd);
    p.paths.data = fullfile(a,[b '_data']);
end

if ~exist(p.paths.data,'dir'); mkdir(p.paths.data); end
data = repmat(struct,length(p.info.subjList),1);
disp(['Getting data' newline 'from: ' p.info.dataURL newline 'to:   ' p.paths.data])
for subjInd = 1:2%length(p.info.subjList)
    if ~exist(fullfile(p.paths.data,[p.info.subjList{subjInd} '.mat']),'file')
        disp(['Downloading subj ' p.info.subjList{subjInd}])
        url = [p.info.dataURL '/files/' p.info.subjList{subjInd} '.mat?download=1'];
        urlwrite(url,fullfile(p.paths.data,[p.info.subjList{subjInd} '.mat']));
    else
        disp([p.info.subjList{subjInd} ' already downloaded'])
    end
end
disp('Downloading done')