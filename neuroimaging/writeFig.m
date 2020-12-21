function writeFig(f,subDir,fileName)
fullDir = fullfile(pwd,subDir);
if ~exist(fullDir,'dir'); mkdir(fullDir); end
fullFileName = fullfile(fullDir,fileName);

f.Color = 'none';
set(findobj(f.Children,'type','Axes'),'color','none')

ext = 'svg';
saveas(f,[fullFileName '.' ext]); disp([fullFileName '.' ext])

f.Color = 'w';

ext = 'fig';
saveas(f,[fullFileName '.' ext]); disp([fullFileName '.' ext])
ext = 'jpg';
saveas(f,[fullFileName '.' ext]); disp([fullFileName '.' ext])