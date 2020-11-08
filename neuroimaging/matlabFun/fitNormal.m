function pd1 = fitNormal(circT,plotIt)
if ~exist('plotIt','var')
    plotIt = 0;
end
%CREATEFIT    Create plot of datasets and fits
%   PD1 = CREATEFIT(CIRCT)
%   Creates a plot, similar to the plot in the main distribution fitting
%   window, using the data that you provide as input.  You can
%   apply this function to the same data you used with dfittool
%   or with different data.  You may want to edit the function to
%   customize the code and this help message.
%
%   Number of datasets:  1
%   Number of fits:  1
%
%   See also FITDIST.

% This function was automatically generated on 12-May-2016 14:46:47

% Output fitted probablility distribution: PD1

% Data from dataset "circT data":
%    Y = circT

% Force all inputs to be column vectors
circT = circT(:);

if plotIt
    % Prepare figure
    clf;
    hold on;
    LegHandles = []; LegText = {};
    
    
    % --- Plot data originally in dataset "circT data"
    [CdfF,CdfX] = ecdf(circT,'Function','cdf');  % compute empirical cdf
    BinInfo.rule = 1;
    [~,BinEdge] = internal.stats.histbins(circT,[],[],BinInfo,CdfF,CdfX);
    [BinHeight,BinCenter] = ecdfhist(CdfF,CdfX,'edges',BinEdge);
    hLine = bar(BinCenter,BinHeight,'hist');
    set(hLine,'FaceColor','none','EdgeColor',[0.333333 0 0.666667],...
        'LineStyle','-', 'LineWidth',1);
    xlabel('Data');
    ylabel('Density')
    LegHandles(end+1) = hLine;
    LegText{end+1} = 'circT data';
    
    % Create grid where function will be computed
    XLim = get(gca,'XLim');
    XLim = XLim + [-1 1] * 0.01 * diff(XLim);
    XGrid = linspace(XLim(1),XLim(2),100);
end

% --- Create fit "fit 1"

% Fit this distribution to get parameter values
% To use parameter estimates from the original fit:
%     pd1 = ProbDistUnivParam('normal',[ 0.132836323, 1.397979560532])
pd1 = fitdist(circT, 'normal');
if plotIt
    YPlot = pdf(pd1,XGrid);
    hLine = plot(XGrid,YPlot,'Color',[1 0 0],...
        'LineStyle','-', 'LineWidth',2,...
        'Marker','none', 'MarkerSize',6);
    LegHandles(end+1) = hLine;
    LegText{end+1} = 'fit 1';
    
    % Adjust figure
    box on;
    hold off;
    
    % Create legend from accumulated handles and labels
    hLegend = legend(LegHandles,LegText,'Orientation', 'vertical', 'Location', 'NorthEast');
    set(hLegend,'Interpreter','none');
end