function pd1 = nonparamDistFit(tmpR,plotFlag)
if ~exist('plotFlag','var')
    plotFlag = 0;
end
%CREATEFIT    Create plot of datasets and fits
%   PD1 = CREATEFIT(TMPR)
%   Creates a plot, similar to the plot in the main distribution fitter
%   window, using the data that you provide as input.  You can
%   apply this function to the same data you used with distributionFitter
%   or with different data.  You may want to edit the function to
%   customize the code and this help message.
%
%   Number of datasets:  1
%   Number of fits:  1
%
%   See also FITDIST.

% This function was automatically generated on 26-Mar-2021 10:53:09

% Output fitted probablility distribution: PD1

% Data from dataset "tmpR data":
%    Y = tmpR

% Force all inputs to be column vectors
tmpR = tmpR(:);

if plotFlag
    % Prepare figure
    figure('WindowStyle','docked');
%     clf;
    hold on;
    LegHandles = []; LegText = {};
end


% --- Plot data originally in dataset "tmpR data"
[CdfF,CdfX] = ecdf(tmpR,'Function','cdf');  % compute empirical cdf
if plotFlag
    
    BinInfo.rule = 1;
    [~,BinEdge] = internal.stats.histbins(tmpR,[],[],BinInfo,CdfF,CdfX);
    [BinHeight,BinCenter] = ecdfhist(CdfF,CdfX,'edges',BinEdge);
    hLine = bar(BinCenter,BinHeight,'hist');
    set(hLine,'FaceColor','none','EdgeColor',[0.333333 0 0.666667],...
        'LineStyle','-', 'LineWidth',1);
    xlabel('Data');
    ylabel('Density')
    LegHandles(end+1) = hLine;
    LegText{end+1} = 'tmpR data';
    
    % Create grid where function will be computed
    XLim = get(gca,'XLim');
    XLim = XLim + [-1 1] * 0.01 * diff(XLim);
    XGrid = linspace(XLim(1),XLim(2),100);
end


% --- Create fit "fit 3"
pd1 = fitdist(tmpR,'kernel','kernel','normal','support','unbounded');
if plotFlag
    YPlot = pdf(pd1,XGrid);
    hLine = plot(XGrid,YPlot,'Color',[0.666667 0.333333 0],...
        'LineStyle','-', 'LineWidth',2,...
        'Marker','none', 'MarkerSize',6);
    LegHandles(end+1) = hLine;
    LegText{end+1} = 'fit 3';
    
    % Adjust figure
    box on;
    hold off;
    
    % Create legend from accumulated handles and labels
    hLegend = legend(LegHandles,LegText,'Orientation', 'vertical', 'FontSize', 9, 'Location', 'northeast');
    set(hLegend,'Interpreter','none');
end
