function out = nii2mat
% [FILENAME_mask, PATHNAME_mask, ~] = uigetfile({'*.nii.gz' '*.nii'}','SelectMaskFile','MultiSelect', 'on');
% [FILENAME_data, PATHNAME_data, ~] = uigetfile({'*.nii.gz' '*.nii'}','SelectDataFile','MultiSelect', 'on');
FILENAME_mask = {'v1.nii.gz'    'v2.nii.gz'    'v3.nii.gz'};
FILENAME_data = {'lh.ecc.nii.gz'    'rh.ecc.nii.gz'};

PATHNAME_mask = uigetdir;
PATHNAME_data = PATHNAME_mask;
% PATHNAME_mask = '/mnt/hgfs/Work/projects/150422_decodingHRF/C_processingComplexRFEfix/01ha/run1c_masking/';
% PATHNAME_data = '/mnt/hgfs/Work/projects/150422_decodingHRF/C_processingComplexRFEfix/01ha/run1c_masking/';

%% Combine masks
for i = 1:length(FILENAME_mask)
    [~, image(:,:,:,i), ~, ~] = BrikLoad(fullfile(PATHNAME_mask,FILENAME_mask{i}));
end
ind = [];
indOver = [];
imageMat = [];
imageMatOver = [];
for z = 1:size(image,3)
    for y = 1:size(image,2)
        for x = 1:size(image,1)
            if any(image(x,y,z,:))
                if length(find(image(x,y,z,:)))==1
                    ind(end+1,:) = [x y z];
                    imageMat(end+1,:) = find(image(x,y,z,:));
                elseif length(find(image(x,y,z,:)))>=1
                    indOver(end+1,:) = [x y z];
                    imageMatOver(end+1,:) = find(image(x,y,z,:));
                end
            end
        end
    end
end

out.mask.im             = imageMat;
out.mask.ind            = ind;
out.mask.path           = fullfile(PATHNAME_mask,FILENAME_mask);
out.mask.overlap.im     = imageMatOver;
out.mask.overlap.ind    = indOver;



mask = image; clear image
%% Combine data
for i = 1:length(FILENAME_data)
    [~, image(:,:,:,i), ~, ~] = BrikLoad(fullfile(PATHNAME_data,FILENAME_data{i}));
end

ind = [];
indOver = [];
imageMat = [];
imageMatOver = [];
for i = 1:size(out.mask.ind,1)
    x = out.mask.ind(i,1); y = out.mask.ind(i,2); z = out.mask.ind(i,3);
    if all(squeeze(image(x,y,z,:))==[0; 0]) || length(find(image(x,y,z,:)))==1
        ind(end+1,:) = [x y z];
        tmpInd = find(image(x,y,z,:));
        if tmpInd
            imageMat(end+1,:) = image(x,y,z,tmpInd);
        else
            imageMat(end+1,:) = 0;
        end
    else
        indOver(end+1,:) = [x y z];
        imageMatOver(end+1,:) = image(x,y,z,:);
    end
end

out.data.im             = imageMat;
out.data.ind            = ind;
out.data.path           = fullfile(PATHNAME_data,FILENAME_data);
out.data.overlap.im     = imageMatOver;
out.data.overlap.ind    = indOver;

%% Save
ecc = out;
save(fullfile(PATHNAME_data,'eccNmask.mat'),'ecc')