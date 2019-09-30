function [fwhm] = fit_resolution(spec,xpk)
% find FWHM of selected peaks in GBS spectra
% INPUT: 
%    spectrum in question
%    position of each peak (in pixels) 
% OUTPUT:
%    FWHM of each peak (in pixels)

% initialize output array
fwhm=zeros(size(xpk));

% loop over all peaks
for i=1:length(xpk)

    ypk=spec(xpk(i));
    
    % assume width is always <22 pixels
    hw=11;

    % left side of peak
    ls=interp1(spec(xpk(i)-hw:xpk(i)),[xpk(i)-hw:xpk(i)],ypk/2);
    % right side of peak
    rs=interp1(spec(xpk(i):xpk(i)+hw),[xpk(i):xpk(i)+hw],ypk/2);

    % FWHM
    fwhm(i)=rs-ls;

end

end
