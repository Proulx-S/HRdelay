clear all
close all

old = load('C:\Users\sebas\OneDrive - McGill University\McGill\work\projects\170210_HRdecoding\C_processing\02jp_sess1.mat');
new = load('C:\Users\sebas\OneDrive - McGill University\dataBig\C-derived\DecodingHR\fun\y\02jp\v1SinCos_1perRun_move12.mat');

size(old.d.normData,2)
numel(new.results.OLS.mixed.delay(:,:,:,1))


hist(old.d.ecc)

new.results.OLS.mixed




%Load ecc mask
                    a = load_nii(fullfile(dataDirIn,'run1c_masking',['lh.ecc.nii.gz']));
                    eccL = flipdim(permute(a.img,[3 1 2 4]),1); clear a
                    eccL(~metric.mask)=0;
                    a = load_nii(fullfile(dataDirIn,'run1c_masking',['rh.ecc.nii.gz']));
                    eccR = flipdim(permute(a.img,[3 1 2 4]),1); clear a
                    eccR(~metric.mask)=0;
                    eccMask = zeros(size(metric.mask));
                    eccMask(eccL>=p.eccLowLim&eccL<p.eccHighLim) = 1;
                    eccMask(eccR>=p.eccLowLim&eccR<p.eccHighLim) = 1;
                    metric.mask(~eccMask) = 0;
