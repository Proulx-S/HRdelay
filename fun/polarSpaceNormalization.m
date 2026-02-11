function [x,polNorm] = polarSpaceNormalization(x,chanSpace,te,phaseOffset)
polNorm.orig = x;
[x,sz] = reDim1(x);
if ~exist('te','var') || isempty(te)
    te = false(size(x,1),1);
end
if isstruct(chanSpace)
    polNorm = chanSpace.pol;
    chanSpace = chanSpace.chanSpace;
    computeNorm = false;
else
    computeNorm = true;
end

switch chanSpace
    case {'cart_affineRot' 'cartImag_affineRot' 'cartReal_affineRot'}
        normSpace.rhoScale = 'vox';
        normSpace.thetaShift = 'roi';
    case 'cartRoi'
        normSpace.rhoScale = 'roi';
        normSpace.thetaShift = 'roi';
    case {'cart' 'cartReal' 'cartImag'}
        normSpace.rhoScale = 'vox';
        normSpace.thetaShift = 'vox';
    case {'cartNoAmp' 'cartNoAmpImag' 'delay'}
        normSpace.rhoScale = 'rm';
        normSpace.thetaShift = 'vox';
    case {'cartNoAmp_affineRot' 'cartNoAmp_affineRot_affineCart' 'cartNoAmpImag_affineRot'}
        normSpace.rhoScale = 'rm';
        normSpace.thetaShift = 'roi';
    case 'cartNoDelay'
        normSpace.rhoScale = 'vox';
        normSpace.thetaShift = 'rm';
    otherwise
        error('X')
end

%% Compute norm
if computeNorm
    % rho scale
    switch normSpace.rhoScale
        case 'vox'
            polNorm.rhoScale = abs(mean(x(~te,:),1));
        case 'roi'
            polNorm.rhoScale = abs(mean(mean(x(~te,:),1),2));
        case 'none'
            polNorm.rhoScale = ones([1 size(x(~te,:),2)]);
        case 'rm'
            polNorm.rhoScale = nan;
        otherwise
            error('X')
    end
    % theta shift
    switch normSpace.thetaShift
        case 'vox'
            polNorm.thetaShift = angle(mean(x(~te,:),1));
            
            if ~strcmp(chanSpace,'delay')
                % identify and preserve neg response voxels
                theta = angle(mean(mean(x(~te,:),1),2));
                indNeg = wrapToPi(polNorm.thetaShift - theta+pi/2)<0;
                polNorm.thetaShift(1,indNeg) = -wrapToPi( pi - polNorm.thetaShift(1,indNeg) );
            end
        case 'roi'
            polNorm.thetaShift = angle(mean(mean(x(~te,:),1),2));
        case 'none'
            polNorm.thetaShift = zeros([1,size(x(~te,:),2)]);
        case 'rm'
            polNorm.thetaShift = nan;
        otherwise
            error('X')
    end
end
if exist('phaseOffset','var')
    polNorm.thetaShift = polNorm.thetaShift + phaseOffset;
end
%% Apply norm
% rho scale
switch normSpace.rhoScale
    case {'vox' 'roi'}
        rho = abs(x)./polNorm.rhoScale;
    case 'none'
        rho = abs(x);
    case 'rm'
        rho = ones(size(x));
    otherwise
        error('X')
end
% theta shift
switch normSpace.thetaShift
    case {'vox' 'roi'}
        theta = wrapToPi(angle(x) - polNorm.thetaShift);
    case 'none'
        theta = angle(x);
    case 'rm'
        theta = zeros(size(x));
    otherwise
        error('X')
end
% theta spread
if ~strcmp(chanSpace,'delay')
    [u,v] = pol2cart(theta,rho); clear theta rho
else
    %tricky trick to keep data packed in complex numbers: put angle in the
    %imaginary part of x
    u = 0;
    v = theta;
end
x = complex(u,v); clear u v

x = reDim2(x,sz);


function [x,sz] = reDim1(x)
sz = size(x);
nDim = length(sz);
switch nDim
    case 2
    case 3
        % Redim
        xTmp = nan(prod(sz([1 3])),sz(2));
        teTmp = nan(prod(sz([1 3])),1);
        for i = 1:sz(3)
            xTmp( (1:sz(1)) + sz(1)*(i-1) , : ) = x(:,:,i);
            teTmp( (1:sz(1)) + sz(1)*(i-1) , 1 ) = false(sz(1),1);
        end
        x = xTmp;
        te = teTmp;
    otherwise
        error('dim1 must be samples and dim2 voxels, that''s it')
end


function x = reDim2(x,sz)
nDim = length(sz);
if nDim==3
    xTmp = nan(sz);
    for i = 1:sz(3)
        xTmp(:,:,i) = x( (1:sz(1)) + sz(1)*(i-1) , : );
    end
    x = xTmp;
end