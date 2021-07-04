function [eccXY,densityXY,X,Y,okInd,U,V,eccXYspline,densityXYspline,eccXYbin,densityXYmed] = getEccPd(d,ecc,sm,plotFlag)
if exist('ecc','var') || ~isempty(ecc)
    d.voxProp.ecc = ecc;
end
if ~exist('plotFlag','var')
    plotFlag = 0;
end
[densityUV,U,V,densityXY,X,Y] = pol2surf(d.voxProp);
if plotFlag>1
    figure('WindowStyle','docked');
    surf(X,Y,densityXY,'LineStyle','none'); hold on
    scatter3(U,V,densityUV,eps,'w.')
    colormap hot
    alpha(.8)
    set(gca, 'color', 'b');
    ax = gca;
    ax.PlotBoxAspectRatio(1:2) = max(ax.PlotBoxAspectRatio(1:2));
    ax.DataAspectRatio(1:2) = max(ax.DataAspectRatio(1:2));
end
[~,eccXY] = cart2pol(X,Y);
okInd = ~isnan(densityXY);
if exist('sm','var') && ~isempty(sm)
    sm2 = sm;
end
sm = 0.8;
%median filter
n = 100;
if ~any(d.voxProp.hemifieldL) || ~any(d.voxProp.hemifieldR)
    n = n/2;
end
[eccXYbin,densityXYmed] = medianFilter(eccXY(okInd),densityXY(okInd),n);

%smothing spline
% [fitresult, ~] = fitSpline(eccXY(okInd), densityXY(okInd), sm);
[fitresult, ~] = fitSpline(eccXYbin, densityXYmed, sm2);
eccXYspline = linspace(min(eccXY(okInd)),max(eccXY(okInd)),100);
densityXYspline = fitresult(eccXYspline);
densityXYspline(densityXYspline<0) = nan;


function [fitresult, gof] = fitSpline(rho, densityXY, sm, plotFlag)
warning('off','curvefit:prepareFittingData:removingNaNAndInf');
if ~exist('plotFlag','var')
    plotFlag = 0;
end

%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( rho, densityXY );

% Set up fittype and options.
ft = fittype( 'smoothingspline' );
opts = fitoptions( 'Method', 'SmoothingSpline' );
if exist('sm','var') && ~isempty(sm)
    opts.SmoothingParam = sm;
else
    opts.SmoothingParam = 0.8;
end

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

if plotFlag
    % Plot fit with data.
    figure( 'Name', 'untitled fit 1' );
    h = plot( fitresult, xData, yData );
    legend( h, 'densityXY vs. rho', 'untitled fit 1', 'Location', 'NorthEast', 'Interpreter', 'none' );
    % Label axes
    xlabel( 'rho', 'Interpreter', 'none' );
    ylabel( 'densityXY', 'Interpreter', 'none' );
    grid on
end
