function compareDataType(res,Xlabel,Ylabel,svmSpace,titleStr)

figure('WindowStyle','docked');
compareRes(res.([svmSpace '_' Xlabel]),res.([svmSpace '_' Ylabel]))
xlabel(Xlabel)
ylabel(Ylabel)
title([svmSpace '; ' titleStr])
