function cal_spec(pixel, wv, out_file_nm, n_pix)
% --------------------------------------------------------------
% Wavelength calibration for spectra
% Cristen Adams
% Written Sep 2010
%
%cal_spec([1999 1929 1782 1660 1536 1508 1310 1129 970], [396.9 393.4 386 379.8 373.5 372.1 361.9 352.5 344.2 ], 'cal_P0_2013_64.clb', 2048)
% Creates a calibration file to convert
% pixels to wavelengths that can be read in WINDOAS based
% on a second order polynomial fit of measured pixel vs wavelength
%
% INPUT:
%   pixel = list of pixels (float or int)
%   wv = list of wavelengths corresponding to pixels
%   n_pix = number of pixels
%   out_file_nm = string output file name
% eg:cal_spec([397,493,526,652,881,924,1383,1441],
% [383.2,393.3,396.8,410.1,434,438.4,486.1,492], 'cal_U1_2012_130.clb', 2048)
% -------------------------------------------------------------

% Check for errors
if length(pixel) ~= length(wv)
    disp('ERROR: Pixcel and wavelength vectors must be the same size.')
end

CCD = 1:n_pix;

% Makes a pixel vs Fraunhofer polynomial fit
C= polyfit(pixel, wv,2);
T1 = [num2str(C(1)) 'X^2 + ' num2str(C(2)) 'X + ' num2str(C(3))];
disp(T1)
calibration = C(1) * CCD .^ 2 + C(2) * CCD + C(3);

% Make a figure
figure
hold on
plot(pixel, wv, 'o')
plot(CCD, calibration)
xlabel('Pixel')
ylabel('Wavelength (nm)')
title(T1)

% Print to file
cal_file = fopen(out_file_nm, 'wt'); % write text mode
fprintf(cal_file, '%f\n', calibration);
fclose(cal_file);