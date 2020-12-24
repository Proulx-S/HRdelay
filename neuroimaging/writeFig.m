function writeFig(f,subDir,fileName,verbose)
if ~exist('verbose','var')
    verbose = 1;
end
fullDir = fullfile(pwd,subDir);
if ~exist(fullDir,'dir'); mkdir(fullDir); end
fullFileName = fullfile(fullDir,fileName);

f.Color = 'none';
set(findobj(f.Children,'type','Axes'),'color','none')

ext = 'svg';
saveas(f,[fullFileName '.' ext]); if verbose; disp([fullFileName '.' ext]); end

f.Color = 'w';

ext = 'fig';
saveas(f,[fullFileName '.' ext]); if verbose; disp([fullFileName '.' ext]); end
ext = 'jpg';
saveas(f,[fullFileName '.' ext]); if verbose; disp([fullFileName '.' ext]); end