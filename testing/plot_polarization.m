% To load polarized spectra and plot results

% path='/home/kristof/work/GBS/UT-GBS/2017/testing_Eureka/';
% path='/home/kristof/work/GBS/PEARL-GBS/2017/testing_Eureka/';
% path='/home/kristof/work/GBS/UT-GBS/2018/testing_Eureka/testing/';
path='/Users/raminaalwarda/Desktop/PhysicsPhD/GBSdata/UT-GBSdata/UTGBS_2020/testing_Eureka/csv/'
% path='./';

% % UTGBS Eureka 2019
% spec_0=get_spectrum_v2([path '28021820.csv'],'8:18:47 PM');
% spec_45=get_spectrum_v2([path '28021820.csv'],'8:21:55 PM');
% spec_90=get_spectrum_v2([path '28021820.csv'],'8:27:30 PM');
% spec_135=get_spectrum_v2([path '28021820.csv'],'8:30:45 PM');
% spec_180=get_spectrum_v2([path '28021820.csv'],'8:33:54 PM');

% % PGBS Eureka 2019
% spec_0=get_spectrum_v2([path '01031919.csv'],'7:10:03 PM');
% spec_45=get_spectrum_v2([path '01031919.csv'],'7:11:12 PM');
% spec_90=get_spectrum_v2([path '01031919.csv'],'7:12:22 PM');
% spec_135=get_spectrum_v2([path '01031919.csv'],'7:13:31 PM');
% spec_180=get_spectrum_v2([path '01031919.csv'],'7:14:40 PM');
% spec_270=get_spectrum_v2([path '01031919.csv'],'7:07:45 PM');
% spec_315=get_spectrum_v2([path '01031919.csv'],'7:08:54 PM');
% 
% PGBS Eureka 2020
% filename1 = '27022016.csv'
% filename2 = '27022017.csv'
% spec_0=get_spectrum_v2([path filename1],'4:53:41 PM');
% spec_45=get_spectrum_v2([path filename1],'4:55:57 PM');
% spec_90=get_spectrum_v2([path filename1],'4:58:13 PM');
% spec_135=get_spectrum_v2([path filename2],'5:00:29 PM');
% spec_180=get_spectrum_v2([path filename2],'5:02:45 PM');
% spec_225=get_spectrum_v2([path filename2],'5:05:01 PM');
% spec_270=get_spectrum_v2([path filename2],'5:07:17 PM');
% % spec_315=get_spectrum_v2([path filename2],'8:40:46 PM');
% % spec_360=get_spectrum_v2([path filename2],'8:43:03 PM');

% UTGBS Eureka 2020
filename1 = '02032018.csv'
filename2 = '02032019.csv'
spec_0=get_spectrum_v2([path filename1],'6:56:00 PM');
spec_45=get_spectrum_v2([path filename2],'7:01:24 PM');
spec_90=get_spectrum_v2([path filename2], '7:05:29 PM');
spec_135=get_spectrum_v2([path filename2],'7:10:58 PM');
spec_180=get_spectrum_v2([path filename2],'7:15:04 PM');
spec_225=get_spectrum_v2([path filename2],'7:19:09 PM');
spec_270=get_spectrum_v2([path filename2],'7:23:15 PM');
spec_315=get_spectrum_v2([path filename2],'7:27:20 PM');
spec_360=get_spectrum_v2([path filename2],'7:31:26 PM');

figure()
plot(spec_0./spec_90, 'r'), hold on
plot(spec_225./spec_315, 'b'), hold on
%plot(spec_0./spec_360, 'g'), hold on
legend('0/90','225/315', location','southeast')
% plot(spec_90./spec_180, 'k'), hold on
xlim([0,2048])
xlabel('CCD Pixels')
ylabel('Ratio of spectra')


% % PGBS Eureka 2020
% filename = '26022020.csv'
% spec_0=get_spectrum_v2([path filename],'8:22:37 PM');
% spec_45=get_spectrum_v2([path filename],'8:26:02 PM');
% spec_90=get_spectrum_v2([path filename],'8:28:18 PM');
% spec_135=get_spectrum_v2([path filename],'8:31:42 PM');
% spec_180=get_spectrum_v2([path filename],'8:33:58 PM');
% spec_225=get_spectrum_v2([path filename],'8:36:14 PM');
% spec_270=get_spectrum_v2([path filename],'8:38:30 PM');
% spec_315=get_spectrum_v2([path filename],'8:40:46 PM');
% spec_360=get_spectrum_v2([path filename],'8:43:03 PM');


% % UTGBS Eureka 2018
% spec_0=get_spectrum_v2([path '28021820.csv'],'8:18:47 PM');
% spec_45=get_spectrum_v2([path '28021820.csv'],'8:21:55 PM');
% spec_90=get_spectrum_v2([path '28021820.csv'],'8:27:30 PM');
% spec_135=get_spectrum_v2([path '28021820.csv'],'8:30:45 PM');
% spec_180=get_spectrum_v2([path '28021820.csv'],'8:33:54 PM');

% % UTGBS Eureka 2018 try 2
% spec_0=get_spectrum_v2([path '01031814.csv'],'2:47:59 PM');
% spec_45=get_spectrum_v2([path '01031814.csv'],'2:54:20 PM');
% spec_90=get_spectrum_v2([path '01031814.csv'],'2:51:17 PM');
% % spec_135=get_spectrum_v2([path '01031814.csv'],'2:57:28 PM');
% spec_135=get_spectrum_v2([path '01031815.csv'],'3:02:40 PM');
% % spec_180=get_spectrum_v2([path '01031814.csv'],'2: PM');


% % PGBS Eureka 2018
% spec_0=get_spectrum_v2([path '26021819.csv'],'7:12:59 PM');
% spec_45=get_spectrum_v2([path '26021819.csv'],'7:16:11 PM');
% spec_90=get_spectrum_v2([path '26021819.csv'],'7:19:01 PM');
% spec_135=get_spectrum_v2([path '26021819.csv'],'7:22:07 PM');
% spec_180=get_spectrum_v2([path '26021819.csv'],'7:24:54 PM');

% UT-GBS, Eureka 2016
% % spec=get_spectrum_v2([path '17031615.csv'],'3:56:05 PM');
% % spec_0=get_spectrum_v2([path '17031616.csv'],'4:01:31 PM');
% % spec_45=get_spectrum_v2([path '17031616.csv'],'4:12:07 PM');
% % spec_90=get_spectrum_v2([path '17031616.csv'],'4:06:53 PM');
% % spec_135=get_spectrum_v2([path '17031616.csv'],'4:18:53 PM');
% % % spec_135=get_spectrum_v2('17031616.csv','4:23:28 PM');
% % % spec_180=get_spectrum_v2('08031618.csv','5:57:52 PM');

% PGBS, Eureka 2017
% % spec=get_spectrum_v2([path '27021717.csv'],'5:06:12 PM');
% % spec_0=get_spectrum_v2([path '27021717.csv'],'5:11:17 PM');
% % spec_45=get_spectrum_v2([path '27021717.csv'],'5:14:33 PM');
% % spec_90=get_spectrum_v2([path '27021717.csv'],'5:17:49 PM');
% % spec_135=get_spectrum_v2([path '27021717.csv'],'5:21:30 PM');
% % spec_180=get_spectrum_v2('27021717.csv','5:25:14 PM');

% % UT-GBS Eureka 2017 (?)
% % spec=get_spectrum_v2([path '27021717.csv'],'5:06:12 PM');
% spec_0=get_spectrum_v2([path '07031720.csv'],'8:40:36 PM');
% spec_45=get_spectrum_v2([path '07031720.csv'],'8:43:33 PM');
% spec_90=get_spectrum_v2([path '07031720.csv'],'8:46:31 PM');
% spec_135=get_spectrum_v2([path '07031720.csv'],'8:49:28 PM');
% % spec_180=get_spectrum_v2('27021717.csv','5:25:14 PM');
