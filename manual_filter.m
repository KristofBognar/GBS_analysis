% filter oubvious outliers that escaped the standard error filters

tg=1;

figure(99)

gscatter(VCD_table.fd,VCD_table.mean_vcd,VCD_table.ampm,'br','..'); hold on

ind=filter_VCD_output( tg, VCD_table, rcd_S );

ind_out=setdiff(1:length(VCD_table.fd),ind)';

gscatter(VCD_table.fd(ind_out),VCD_table.mean_vcd(ind_out),...
         VCD_table.ampm(ind_out),'br','oo',9); hold on

legend('Unfiltered am','Unfiltered pm','Removed am','Removed pm')


% UT-GBS no2
% filt_1999=[28]; too high, other twilight filtered as well
% filt_2000=[];
% filt_2003=[13:end]; based on Cristen's lack of data there
% filt_2004=[];
% filt_2005=[];
% filt_2006=[2,45]; too high, other twilight filtered as well
% filt_2007=[38]; too low, filtered in Cristen's file
% filt_2008=[];
% filt_2009=[3,24,26]; surrounded by bad values
% filt_2010=[174]; too low, other twilight filtered as well
% filt_2011=[];
% filt_2012=[134]; too low, other twilight filtered as well
% filt_2013=[180]; too low, other twilight filtered as well
% filt_2014=[372]; too high, other twilight filtered as well
% filt_2015=[195]; too few twilight spectra
% filt_2016=[];
% filt_2017=[];

% UT-GBS o3
% filt_1999=[];
% filt_2000=[];
% filt_2003=[];
% filt_2004=[];
% filt_2005=[17]; too high, other twilight filtered as well
% filt_2006=[];
% filt_2007=[27]; other twilight missing
% filt_2008=[];
% filt_2009=[7]; too low, other twilight filtered as well
% filt_2010=[91]; too high
% filt_2011=[2,453]; outliers, other twilight filtered as well
% filt_2012=[];
% filt_2013=[178]; too low, other twilight missing
% filt_2014=[154,404,405]; too high, other twilight filtered as well; sun too low
% filt_2015=[3]; other twilight filteredas well; 
% filt_2016=[9]; other twilight missing
% filt_2017=[];
