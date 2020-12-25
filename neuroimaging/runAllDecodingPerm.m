function res = runAllDecodingPerm(res,figOption,verbose)
%% Permutation tests
nPerm = 2^14; % will not run if already run, and instead just load previous results

svmSpace = 'cart';
res.(svmSpace) = runDecoding(res.(svmSpace),verbose,nPerm);

svmSpace = 'cartNoAmp';
res.(svmSpace) = runDecoding(res.(svmSpace),verbose,nPerm);

svmSpace = 'cartNoDelay';
res.(svmSpace) = runDecoding(res.(svmSpace),verbose,nPerm);

svmSpace = 'cartReal';
res.(svmSpace) = runDecoding(res.(svmSpace),verbose,nPerm);

%% Plot everything
accPerm = plotDecodingPerm_acc(res,figOption,verbose);
plotDecoding_acc(res,figOption,verbose,accPerm)
