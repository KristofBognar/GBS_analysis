function plot_straylight()
% Read and plot stray light compared to regular spectrum 
% based on straylight tests
%
% INPUT: Prompt user for spectra files and start times
% OUTPUT: straylight / regular spectrum plot


% specfile='06091613.csv';
% spectime='1:26:17';
% dcfile='06091613.csv';
% dctime='1:24:42';
% redfile='06091613.csv';
% redtime='1:21:33';

% specfile='06091613.csv';
% spectime='1:08:38';
% dcfile='06091613.csv';
% dctime='1:10:03';
% redfile='06091613.csv';
% redtime='1:12:52';

% PGBS Eureka 2017 run 1
% specfile='27021717.csv';
% spectime='5:33:30';
% dcfile=specfile;
% dctime='5:57:28';
% redfile=specfile;
% redtime='5:37:29';

% % PGBS Eureka 2017 run 2
% specfile='28021716.csv';
% spectime='4:22:19';
% dcfile='28021715.csv';
% dctime='3:38:19';
% redfile=specfile;
% redtime='4:25:59';
% 
% % UT-GBS Eureka 2017 
% specfile='07031720.csv';
% spectime='8:20:35';
% dcfile=specfile;
% dctime='8:22:04';
% redfile=specfile;
% redtime='8:24:57';

% % PGBS Eureka 2018 run 4
% specfile='27021815.csv';
% spectime='3:54:46';
% dcfile='26021820.csv';
% dctime='8:15:59';
% redfile='27021816.csv';
% redtime='4:02:22';

% % UTGBS Eureka 2018
% specfile='28021820.csv';
% spectime='8:04:35';
% dcfile=specfile;
% dctime='8:06:30';
% redfile=specfile;
% redtime='8:10:02';

% % UTGBS Eureka 2018 try 2
% specfile='01031814.csv';
% spectime='2:33:06';
% dcfile=specfile;
% dctime='2:35:02';
% redfile=specfile;
% redtime='2:38:35';

% PGBS Eureka 2019
% specfile='02031916.csv';
% spectime='4:45:49';
% dcfile='01031921.csv';
% dctime='9:12:17';
% redfile=specfile;
% redtime='4:54:30';

% PGBS Eureka 2020
specfile='27022019.csv';
spectime='7:43:11';
dcfile='27022019.csv';
dctime='7:49:59';
redfile=specfile;
redtime='7:45:27';

% UTGBS Eureka 2020
specfile='02032018.csv';
spectime='6:26:21';
dcfile='02032020.csv';
dctime='8:09:44';
redfile=specfile;
redtime='6:34:26';


path='/Users/raminaalwarda/Desktop/PhysicsPhD/GBSdata/UT-GBSdata/UTGBS_2020/testing_Eureka/csv/';

spec=get_spectrum_v2([path specfile],spectime);
dc=get_spectrum_v2([path dcfile],dctime);

red=get_spectrum_v2([path redfile],redtime);

% load straylight.mat

straylight=red-dc;

plot(straylight./(spec))
xlim([0,2048])
xlabel('CCD Pixel')
ylabel('Straylight / regular spectrum')

ylim([0,0.29])

