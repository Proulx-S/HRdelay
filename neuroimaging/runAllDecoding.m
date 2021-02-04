function res = runAllDecoding(figOption,verbose)
if ~exist('verbose','var')
    verbose = 1;
end

if ~exist('figOption','var') || isempty(figOption)
    figOption.save = 0;
    figOption.subj = 1; % 'all' or subjInd
end

svmSpace = 'cart';
dataType = 'wave';
[resBS.([svmSpace '_' dataType]),resWS.([svmSpace '_' dataType])] = runDecoding(svmSpace,dataType,verbose);
dataType = 'sin';
[resBS.([svmSpace '_' dataType]),resWS.([svmSpace '_' dataType])] = runDecoding(svmSpace,dataType,verbose);
compareDataType(resBS,'wave','sin',svmSpace,'between-session')
compareDataType(resWS,'wave','sin',svmSpace,'within-session')


svmSpace = 'cartNoAmp';
dataType = 'wave';
[resBS.([svmSpace '_' dataType]),resWS.([svmSpace '_' dataType])] = runDecoding(svmSpace,dataType,verbose);
dataType = 'sin';
[resBS.([svmSpace '_' dataType]),resWS.([svmSpace '_' dataType])] = runDecoding(svmSpace,dataType,verbose);
compareDataType(resBS,'wave','sin',svmSpace,'between-session')
compareDataType(resWS,'wave','sin',svmSpace,'within-session')

svmSpace = 'cartNoDelay';
dataType = 'wave';
[resBS.([svmSpace '_' dataType]),resWS.([svmSpace '_' dataType])] = runDecoding(svmSpace,dataType,verbose);
dataType = 'sin';
[resBS.([svmSpace '_' dataType]),resWS.([svmSpace '_' dataType])] = runDecoding(svmSpace,dataType,verbose);
compareDataType(resBS,'wave','sin',svmSpace,'between-session')
compareDataType(resWS,'wave','sin',svmSpace,'within-session')

svmSpace = 'cartReal';
dataType = 'wave';
[resBS.([svmSpace '_' dataType]),resWS.([svmSpace '_' dataType])] = runDecoding(svmSpace,dataType,verbose);
dataType = 'sin';
[resBS.([svmSpace '_' dataType]),resWS.([svmSpace '_' dataType])] = runDecoding(svmSpace,dataType,verbose);
compareDataType(resBS,'wave','sin',svmSpace,'between-session')
compareDataType(resWS,'wave','sin',svmSpace,'within-session')

% [P,~,STATS] = signrank(mean(res.cart.acc,2),mean(res.cartReal.acc,2),'tail','right');
% disp('cart VS cartReal:')
% disp(['signed rank = ' num2str(STATS.signedrank)])
% disp(['one-sided p = ' num2str(P)])

disp('Between-session')
tmpField = fields(resBS);
tmp = struct2cell(resBS);
tmp = cell2struct(tmp(1:2:end),tmpField(1:2:end));
plotDecoding_acc(tmp,figOption,verbose)
disp('Witheen-session')
tmpField = fields(resWS);
tmp = struct2cell(resWS);
tmp = cell2struct(tmp(1:2:end),tmpField(1:2:end));
plotDecoding_acc(tmp,figOption,verbose)

res.BS = resBS;
res.WS = resWS;


% resBS.cartNoAmp.sess.acc
% 1-resBS.cartNoAmp.sess.acc_p