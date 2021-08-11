function diaryON(p)
if p.termOption.save
    outFile = fullfile(p.wd,['termOutputs-' datestr(now,'YYYYmmDD-hh_MM_ss') '.log''']);
    cmd = ['diary ''' outFile];
    eval(cmd)
    disp(datestr(now))
end