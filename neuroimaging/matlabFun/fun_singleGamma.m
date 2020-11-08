function gamma = fun_singleGamma(t,HRF_PARAMETERS)
%Worsley

% HRF_PARAMETERS is a matrix whose rows are 5 parameters for the 
% hemodynamic response function, one row for each event type and column
% of S (if there is just one row, this is repeated as necessary). 
% The hrf is modeled as the difference of two gamma density functions 
% (Glover, NeuroImage, 9:416-429). The components of HRF_PARAMETERS are:
% 1. PEAK1: time to the peak of the first gamma density;
% 2. FWHM1: approximate FWHM of the first gamma density;
% 3. PEAK2: time to the peak of the second gamma density;
% 4. FWHM2: approximate FWHM of the second gamma density;
% 5. DIP: coefficient of the second gamma density;
%    Final hrf is:   gamma1/max(gamma1)-DIP*gamma2/max(gamma2)
%    scaled so that its total integral is 1. 
% If PEAK1=0 then there is no smoothing of that event type with the hrf.
% If PEAK1>0 but FWHM1=0 then the design is simply lagged by PEAK1.
% Default is: [5.4 5.2 10.8 7.35 0.35] chosen by Glover (1999) for 
% an auditory stimulus. 
% If HRF_PARAMETERS is a structure, then HRF_PARAMETERS.T is a matrix
% whose rows are the times in seconds of a user-supplied HRF, one 
% row for each event type and column of S (if there is just one row, 
% this is repeated as necessary).  Times must start at 0 but need not 
% be equally spaced; spacing of 0.02s is recommended. HRF_PARAMETERS.H 
% are the corresponding values of the HRF at those times. 

if ~exist('HRF_PARAMETERS','var') || isempty(HRF_PARAMETERS)
    HRF_PARAMETERS = [5.4 5.2 10.8 7.35 0.35];
elseif length(HRF_PARAMETERS)==1 && HRF_PARAMETERS==1
    HRF_PARAMETERS = [5.4 5.2 0 0 0];
end

peak1=HRF_PARAMETERS(1);
fwhm1=HRF_PARAMETERS(2);
peak2=HRF_PARAMETERS(3);
fwhm2=HRF_PARAMETERS(4);
dip  =HRF_PARAMETERS(5);


alpha=peak1^2/fwhm1^2*8*log(2);
beta=fwhm1^2/peak1/8/log(2);
gamma1=(t/peak1).^alpha.*exp(-(t-peak1)./beta);
if any(t<0)
    gamma1(t<0)=0;
end

if peak2
    alpha=peak2^2/fwhm2^2*8*log(2);
    beta=fwhm2^2/peak2/8/log(2);
    gamma2=(t/peak2).^alpha.*exp(-(t-peak2)./beta);
    if any(t<0)
        gamma2(t<0)=0;
    end
    gamma = gamma1-dip.*gamma2;
else
    gamma = gamma1;
end

gamma = gamma./max(gamma);

