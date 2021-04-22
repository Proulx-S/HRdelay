function [featSel_areaAndFov,d,p] = empiricalFov(d,p,outPath)
filePath = fullfile(outPath,mfilename);
forceFlag = 0;

%% Initiate params
p.featSel.fov.empirical.padFac             = 1.2;
p.featSel.fov.empirical.minContPercentArea = 0.05;
p.featSel.fov.empirical.auto(1).smList           = 0.001; % ecc
p.featSel.fov.empirical.auto(1).mergeRadiusList  = 0.70; % ecc
p.featSel.fov.empirical.auto(1).marginRadiusList = 0.40; % ecc
p.featSel.fov.empirical.auto(2).smList           = 0.25; % ecc
p.featSel.fov.empirical.auto(2).mergeRadiusList  = 0.70; % ecc
p.featSel.fov.empirical.auto(2).marginRadiusList = 0.40; % ecc

if ~exist([filePath '.mat'],'file') || forceFlag
    %% Precompute flattened voxel ecc distribution on fov and delay map
    disp('Flattening ecc dist: computing hemiL')
    voxProp.L = flattenEccDist(d,'L',p,1);
    disp('Flattening ecc dist: computing hemiR')
    voxProp.R = flattenEccDist(d,'R',p,1);
    disp('Flattening ecc dist: done')
    
    
    disp('Delay map for contour: computing hemiL')
    cont.L = prepareDelayFovContour(d,voxProp.L,p);
    disp('Delay map for contour: computing hemiR')
    cont.R = prepareDelayFovContour(d,voxProp.R,p);
    disp('Delay map for contour: done')
    
    %reorganize voxProp
    voxProp2 = cell(size(d));
    for subjInd = 1:size(d,1)
        for sessInd = 1:size(d,2)
            voxProp2{subjInd,sessInd}.L = voxProp.L{subjInd};
            voxProp2{subjInd,sessInd}.R = voxProp.R{subjInd};
        end
        voxProp.L{subjInd} = {};
        voxProp.R{subjInd} = {};
    end
    voxProp = voxProp2; clear voxProp2
    
    
    %% Get empirical fov
    featSel_areaAndFov = cell(size(d));
    f.L = cell(size(d));
    f.R = cell(size(d));
    % disp('computing feature selection stats')
    for subjInd = 1:size(d,1)
        disp(['subj:' num2str(subjInd) '/' num2str(size(d,1))])
        for sessInd = 1:size(d,2)
            p.subjInd = subjInd;
            p.sessInd = sessInd;
            hemi = 'L';
            [featValLR.(hemi),featMethodLR.(hemi),featIndInLR.(hemi),featInfoLR.(hemi),f.(hemi){subjInd,sessInd}] = getFeatSel_areaAndFov(cont.(hemi){subjInd,sessInd},voxProp{subjInd,sessInd}.(hemi),p);
            hemi = 'R';
            [featValLR.(hemi),featMethodLR.(hemi),featIndInLR.(hemi),featInfoLR.(hemi),f.(hemi){subjInd,sessInd}] = getFeatSel_areaAndFov(cont.(hemi){subjInd,sessInd},voxProp{subjInd,sessInd}.(hemi),p);
            
            % Combine hemifields
            featVal = nan(size(voxProp{subjInd,sessInd}.(hemi).hemifield));
            featIndIn = nan(size(voxProp{subjInd,sessInd}.(hemi).hemifield));
            hemi = 'L';
            featVal(voxProp{subjInd,sessInd}.(hemi).hemifield) = featValLR.(hemi);
            featIndIn(voxProp{subjInd,sessInd}.(hemi).hemifield) = featIndInLR.(hemi);
            hemi = 'R';
            featVal(voxProp{subjInd,sessInd}.(hemi).hemifield) = featValLR.(hemi);
            featIndIn(voxProp{subjInd,sessInd}.(hemi).hemifield) = featIndInLR.(hemi);
            featMethod = featMethodLR.L;
            featInfo = featInfoLR.L;
            
            featSel_areaAndFov{subjInd,sessInd}.featVal = featVal; clear featVal featValLR;
            featSel_areaAndFov{subjInd,sessInd}.featIndIn = featIndIn; clear featIndIn featIndInLR;
            featSel_areaAndFov{subjInd,sessInd}.featMethod = featMethod; clear featMethod featMethodLR;
            featSel_areaAndFov{subjInd,sessInd}.featInfo = featInfo; clear featInfo featInfoLR;
        end
    end
    
    fIndList = [1 2 4 7 9 10];
    supTitleList = {'removing small islands' '1st contours' '1st contours processing' '2nd contours' '2nd contours processing' 'Final contours'};
    % fAll.L = cell(size(fIndList));
    % fAll.R = cell(size(fIndList));
    fAll = cell(size(fIndList));
    % supTitleList = {'Contour Definition' 'Contour Processing' 'Final Contour' 'Contour Masking'};
    for i = 1:length(fIndList)
        fInd = fIndList(i);
        %     fAll.L{i} = figure('WindowStyle','docked');
        fAll{i} = figure('WindowStyle','docked');
        [ha, pos] = tight_subplot(size(d,2), size(d,1)*2, 0, 0.1, 0); delete(ha);
        for subjInd = 1:size(d,1)
            for sessInd = 1:size(d,2)
                hemi = 'L';
                ax.(hemi) = copyobj(f.(hemi){subjInd,sessInd}(fInd).Children,fAll{i});
                ax.(hemi).DataAspectRatioMode = 'auto';
                ax.(hemi).PlotBoxAspectRatioMode = 'auto';
                %           ax = copyobj(f.(hemi){subjInd,sessInd}(fInd).Children,fAll.(hemi){i});
                %             ax.Position = pos{(sessInd-1)*size(d,1)+subjInd};
                ax.(hemi).Position = pos{(sessInd-1)*(size(d,1)*2)+(subjInd*2-1)};
                ax.(hemi).Colormap = f.(hemi){subjInd,sessInd}(fInd).Children.Colormap;
                drawnow
                
                hemi = 'R';
                ax.(hemi) = copyobj(f.(hemi){subjInd,sessInd}(fInd).Children,fAll{i});
                ax.(hemi).DataAspectRatioMode = 'auto';
                ax.(hemi).PlotBoxAspectRatioMode = 'auto';
                ax.(hemi).Position = pos{(sessInd-1)*(size(d,1)*2)+(subjInd*2-1)+1};
                ax.(hemi).Colormap = f.(hemi){subjInd,sessInd}(fInd).Children.Colormap;
                drawnow
                
                yLim = [-1 1].*max(abs([ax.L.YLim ax.R.YLim ax.L.XLim(1) ax.R.XLim(2)]));
                xLim = yLim(2);
                ax.L.YLim = yLim;
                ax.R.YLim = yLim;
                ax.L.XLim = [-xLim 0];
                ax.R.XLim = [0 xLim];
                drawnow
                
                ax.L.PlotBoxAspectRatio = [0.5 1 1];
                ax.R.PlotBoxAspectRatio = [0.5 1 1];
                
                ax.L.YAxis.Visible = 'off';
                ax.R.YAxis.Visible = 'off';
            end
        end
        suptitle(supTitleList{i})
        
        %     fInd = fIndList(i);
        %     fAll.R{i} = figure('WindowStyle','docked');
        %     [ha, pos] = tight_subplot(size(d,2), size(d,1), 0, 0.1, 0); delete(ha);
        %     for subjInd = 1:size(d,1)
        %         for sessInd = 1:size(d,2)
        %             ax = copyobj(f.R{subjInd,sessInd}(fInd).Children,fAll.R{i});
        %             ax.Position = pos{(sessInd-1)*size(d,1)+subjInd};
        %             ax.Colormap = f.R{subjInd,sessInd}(fInd).Children.Colormap;
        %             drawnow
        % %             delete(f{subjInd,sessInd}(fInd).Children);
        %         end
        %     end
        %     suptitle(supTitleList{i})
    end
    f = fAll{end};
    save(fullfile(outPath,mfilename),'featSel_areaAndFov','voxProp','cont','fAll','f')
else
    load(fullfile(outPath,mfilename),'featSel_areaAndFov','cont','voxProp')
    if p.figOption.verbose>=1
        load(fullfile(outPath,mfilename),'f')    
    end
end

%% Pack cont into featSel_areaAndFov and voxProp into d
for subjInd = 1:size(d,1)
    for sessInd = 1:size(d,2)
        featSel_areaAndFov{subjInd,sessInd}.cont.L = cont.L{subjInd,sessInd};
        featSel_areaAndFov{subjInd,sessInd}.cont.R = cont.R{subjInd,sessInd};
        d{subjInd,sessInd}.voxProp.L = voxProp{subjInd,sessInd}.L;
        d{subjInd,sessInd}.voxProp.R = voxProp{subjInd,sessInd}.R;

    end
end

