function [x,cartNorm] = cartSpaceNormalization(x,SVMspace,te,doCartSpaceScale)
[x,sz] = reDim1(x);
if ~exist('te','var') || isempty(te)
    te = false(size(x,1),1);
end
if isstruct(SVMspace)
    cartNorm = SVMspace.svm;
    SVMspace = SVMspace.svmSpace;
    computeNorm = false;
else
    computeNorm = true;
end

if computeNorm
    % Specify norm
    switch SVMspace
        case {'cart' 'cartNoAmp'}
            if doCartSpaceScale
                normSpace.realShift = 'vox';
                normSpace.imagShift = 'vox';
                normSpace.realScale = 'vox';
                normSpace.imagScale = 'vox';
            else
                normSpace.realShift = 'roi';
                normSpace.imagShift = 'roi';
                normSpace.realScale = 'none';
                normSpace.imagScale = 'none';
            end
        case {'cartNoDelay' 'cartReal'}
            normSpace.realShift = 'vox';
            normSpace.imagShift = 'rm';
            if doCartSpaceScale
                normSpace.realShift = 'vox';
                normSpace.imagShift = 'rm';
                normSpace.realScale = 'vox';
                normSpace.imagScale = 'rm';
            else
                normSpace.realShift = 'roi';
                normSpace.imagShift = 'roi';
                normSpace.realScale = 'none';
                normSpace.imagScale = 'rm';
            end
        case {'cartNoAmpImag' 'cartImag'}
            if doCartSpaceScale
                normSpace.realShift = 'rm';
                normSpace.imagShift = 'vox';
                normSpace.realScale = 'rm';
                normSpace.imagScale = 'vox';
            else
                normSpace.realShift = 'roi';
                normSpace.imagShift = 'roi';
                normSpace.realScale = 'rm';
                normSpace.imagScale = 'none';
            end
        otherwise
            error('X')
    end
    
    % Compute norm
    u = real(x(~te,:));
    v = imag(x(~te,:));
    %real shift
    switch normSpace.realShift
        case 'roi'
            cartNorm.realShift = mean(u(:));
        case 'vox'
            cartNorm.realShift = mean(u,1);
        case 'none'
            cartNorm.realShift = 0;
        case 'rm'
            cartNorm.realShift = nan;
        otherwise
            error('X')
    end
    %imag shift
    switch normSpace.imagShift
        case 'roi'
            cartNorm.imagShift = mean(v(:));
        case 'vox'
            cartNorm.imagShift = mean(v,1);
        case 'none'
            cartNorm.imagShift = 0;
        case 'rm'
            cartNorm.imagShift = nan;
        otherwise
            error('X')
    end
    %real scale
    switch normSpace.realScale
        case 'roi'
            cartNorm.realScale = std(mean(u,2),[],1);
        case 'vox'
            cartNorm.realScale = std(u,[],1);
        case 'none'
            cartNorm.realScale = 1;
        case 'rm'
            cartNorm.realScale = nan;
        otherwise
            error('X')
    end
    %imag scale
    switch normSpace.imagScale
        case 'roi'
            cartNorm.imagScale = std(mean(v,2),[],1);
        case 'vox'
            cartNorm.imagScale = std(v,[],1);
        case 'none'
            cartNorm.imagScale = 1;
        case 'rm'
            cartNorm.imagScale = nan;
        otherwise
            error('X')
    end
end

% Apply norm
u = real(x);
v = imag(x);
%real shift
if ~isnan(cartNorm.realShift)
    u = u - cartNorm.realShift;
else
    u = 0;
end
% switch normSpace.realShift
%     case {'roi' 'vox' 'none'}
%         u = u - cartNorm.realShift;
%     case 'rm'
%         u = 0;
%     otherwise
%         error('X')
% end
%imag shift
if ~isnan(cartNorm.imagShift)
    v = v - cartNorm.imagShift;
else
    v = 0;
end
% switch normSpace.imagShift
%     case {'roi' 'vox' 'none'}
%         v = v - cartNorm.imagShift;
%     case 'rm'
%         v = 0;
%     otherwise
%         error('X')
% end
%real scale
if ~isnan(cartNorm.realScale)
    u = u ./ cartNorm.realScale;
else
    u = 0;
end
% switch normSpace.realScale
%     case {'roi' 'vox' 'none'}
%         u = u ./ cartNorm.realScale;
%     case 'rm'
%         u = 0;
%     otherwise
%         error('X')
% end
%imag scale
if ~isnan(cartNorm.imagScale)
    v = v ./ cartNorm.imagScale;
else
    v = 0;
end
% switch normSpace.imagScale
%     case {'roi' 'vox' 'none'}
%         v = v ./ cartNorm.imagScale;
%     case 'rm'
%         v = 0;
%     otherwise
%         error('X')
% end
x = complex(u,v);

x = reDim2(x,sz);