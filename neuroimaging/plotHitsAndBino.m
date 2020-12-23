function plotHitsAndBino(hit,n)
% --- Plot data originally in dataset "hit data"
[CdfF,CdfX] = ecdf(hit,'Function','cdf');  % compute empirical cdf
BinEdge = (0:n+1)-0.5;
[BinHeight,BinCenter] = ecdfhist(CdfF,CdfX,'edges',BinEdge);
hLine = bar(BinCenter,BinHeight,'hist'); hold on
hLine.FaceColor = 'k';
xlim([0 n])
xlabel('Hit count');
ylabel('Density')

% Create grid where function will be computed
XLim = get(gca,'XLim');
XLim = XLim + [-1 1] * 0.01 * diff(XLim);
XIncr = max(1,floor((XLim(2)-XLim(1))/100));
XGrid = floor(XLim(1)):XIncr:ceil(XLim(2));


% --- Create fit "fit 1"

% Fit this distribution to get parameter values
% To use parameter estimates from the original fit:
%     pd1 = ProbDistUnivParam('binomial',[ 142, 0.5])
YPlot = binopdf(XGrid,n,0.5);
% pd1 = fitdist(hit, 'binomial', 'n', n);
% YPlot = pdf(pd1,XGrid);
plot(XGrid,YPlot,'Color','r');