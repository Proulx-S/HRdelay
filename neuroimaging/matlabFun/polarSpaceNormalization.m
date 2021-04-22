function [x,polNorm] = polarSpaceNormalization(x,SVMspace,te,phaseOffset)
[x,sz] = reDim1(x);
if ~exist('te','var') || isempty(te)
    te = false(size(x,1),1);
end
if isstruct(SVMspace)
    polNorm = SVMspace.pol;
    SVMspace = SVMspace.svmSpace;
    computeNorm = false;
else
    computeNorm = true;
end

switch SVMspace
    case {'cart_affineRot' 'cartImag_affineRot' 'cartReal_affineRot'}
        normSpace.rhoScale = 'vox';
        normSpace.thetaShift = 'roi';
    case 'cart_roi'
        normSpace.rhoScale = 'roi';
        normSpace.thetaShift = 'roi';
    case {'cart' 'cartReal' 'cartImag'}
        normSpace.rhoScale = 'vox';
        normSpace.thetaShift = 'vox';
    case {'cartNoAmp' 'cartNoAmpImag'}
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
[u,v] = pol2cart(theta,rho); clear theta rho
x = complex(u,v); clear u v

x = reDim2(x,sz);




