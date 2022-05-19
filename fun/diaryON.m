function diaryON(p)
outFile = fullfile(p.termOption.outDir,['termOutputs-' datestr(now,'YYYYmmDD-hh_MM_ss') '.log''']);
cmd = ['diary ''' outFile];
eval(cmd)
disp(datestr(now))
