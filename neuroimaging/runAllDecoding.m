function runAllDecoding(figOption,verbose)
if ~exist('verbose','var')
    verbose = 1;
end

if ~exist('figOption','var') || isempty(figOption)
    figOption.save = 0;
    figOption.subj = 1; % 'all' or subjInd
end

svmSpace = 'cart';
resTmp = runDecoding(svmSpace,verbose);
res.(svmSpace) = resTmp;

svmSpace = 'cartNoAmp';
resTmp = runDecoding(svmSpace,verbose);
res.(svmSpace) = resTmp;

svmSpace = 'cartNoDelay';
resTmp = runDecoding(svmSpace,verbose);
res.(svmSpace) = resTmp;

svmSpace = 'cartReal';
resTmp = runDecoding(svmSpace,verbose);
res.(svmSpace) = resTmp;

[P,H,STATS] = signrank(mean(res.cart.acc,2),mean(res.cartReal.acc,2),'tail','right');
disp('cart VS cartReal:')
disp(['signed rank = ' num2str(STATS.signedrank)])
disp(['one-sided p = ' num2str(P)])

% plotDecoding_acc(res,figOption,verbose)
% return

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
