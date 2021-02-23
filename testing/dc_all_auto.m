function dc_all_auto(filename_orig,biasfile,dcfile,path,flag_plot)

% this function calculates the bias and dark current from a series of
% files, beginning with filename_orig.
%
% call: dc_all_auto(filename,biasfile,dcfile,flag_plot)
% Xiaoyi's note: dc_all_auto('29021200.csv')

% inputs:   filename_orig  - the file that contains the 30 s exposure time
%           biasfile - the output file for the bias.  if no name is
%                      specified, this will be bias.inp.
%           dcfile - the output file for the dark current.  if no name is
%                    specified, this will be dark_c.inp.
%           flag_plot - if set to 1, will plot the dark signal vs. exposure
%                       time for pixel 1000.  if set to 0, will not plot.
%                       if no value is specified, the default is 0.
% outputs:  bias - the bias values
%           dark_current - the dark current values
%           mean_T_all - the mean CCD temperature for all the runs.
%
% Annemarie Fraser, September 2004.
% Modified Cristen Adams, August 2010

num_pixels=2048;     % the number of pixels across the CCD.


% specifies the defaults if there are input arguments missing.
if nargin==1
    biasfile='bias.inp';
    dcfile='dark_c.inp';
    flag_plot=0;
elseif nargin==2
    dcfile='dark_c.inp';
    flag_plot=0;
elseif nargin==3
    flag_plot=0;
end

filename=filename_orig;
dark_sp = [];
Tmean =[];
exp_time = [];
for i=1:7   % loop for going through each file
   
    % this loop finds the next file.
    test=str2num(filename(7:8));
    if i==1  && exist(filename_orig)==2
        filename=filename_orig;
    elseif test<9
        filename_test=[filename(1:7) int2str(str2num(filename(8))+1) filename(9:12)];
        if exist(filename_test)==0 
            filename=[filename(1:7) int2str(str2num(filename(8))+2) filename(9:12)];
        else
            filename=filename_test;
        end
    elseif test>=9 && test<23
        filename_test=[filename(1:6) int2str(str2num(filename(7:8))+1) filename(9:12)];
        if exist(filename_test)==0 
            filename=[filename(1:6) int2str(str2num(filename(7:8))+2) filename(9:12)];
        else
            filename=filename_test;
        end
    elseif test==23
        date_vec_prev = [str2num(['20' filename(5:6)]),...
            str2num(filename(3:4)), str2num(filename(1:2))];
        date_vec_next = datevec(datenum(date_vec_prev) + 1);
        
        yy_str = int2str(date_vec_next(1));
        mm = date_vec_next(2);
        dd = date_vec_next(3);
        if dd < 10
            dd_str = ['0' int2str(dd)];
        else
            dd_str = int2str(dd);
        end
        if mm < 10
            mm_str = ['0' int2str(mm)];
        else
            mm_str = int2str(mm);
        end
        
        filename=[dd_str mm_str yy_str(3:4) '00' filename(9:12)];
    end
    if exist(filename)==0, continue, end

    disp(['Reading file: ' filename])
    

    % this loop goes through the file(i) and creates a matrix of all the intensities in the file.
    for k=1:20
        [data]=textread([path, filename],'',num_pixels,'delimiter',',','headerlines',3+(num_pixels+1)*(k-1));
        if isempty(data),break,end
        [start_date start_time end_date end_time shutter ideal_counts slit grating turret blaze target exp_time_i accumulations ...
                mean_T min_T max_T mean_Tbox gps_signal satellites altitude longitude latitude]=textread([path filename],...
                '%q %q %q %q %q %q %q %q %q %q %q %f %q %f %q %q %q %q %q %q %q %q',1,'delimiter',',','headerlines',2+(num_pixels+1)*(k-1));
        if str2num(cell2mat(shutter)) == 0
            exp_time = [exp_time; exp_time_i];
            Tmean=[Tmean; mean_T];
            dark_sp=[dark_sp data];
        end
    end
end

% Now rearrange so that we only have values after high exposure time
[max_val, max_ind] = max(exp_time);
[min_val, min_ind] = min(exp_time(max_ind:length(exp_time)));
if isempty(max_val),
    disp('ERROR: No DC measurements found in files')
    return
end
if max_val == min_val
    disp('ERROR: No range in exposure times - likely no DC test taken in time range!!!')
    return
end
min_ind = min_ind + max_ind -1;

dark_sp = dark_sp(:, max_ind:min_ind);
exp_time = exp_time(max_ind:min_ind);
Tmean = Tmean(max_ind:min_ind);


% finds the bias and dark current.
for i = 1:num_pixels
    temp=polyfit(exp_time',dark_sp(i,:),1);
    dc(i) = temp(1)'*1000;
    bias(i) = temp(2)';
end
bias = bias';
dc = dc';

figure
subplot(3,1,1)
hold on
plot(exp_time,dark_sp(1,:),'+','color','b')
plot(exp_time,dark_sp(1000,:),'+','color','r')
plot(exp_time,dark_sp(1500,:),'+','color','g')
legend('Pixel 1', 'Pixel 1000', 'Pixel 1500')
xlabel('Exposure time (s)')
ylabel('Dark spectrum count')
plot(exp_time,bias(1)+dc(1)*exp_time/1000,'b','HandleVisibility','Off');
plot(exp_time,bias(1000)+dc(1000)*exp_time/1000,'r','HandleVisibility','Off');
plot(exp_time,bias(1500)+dc(1500)*exp_time/1000,'g','HandleVisibility','Off');
subplot(3,1,2)
plot(bias)
xlabel('Pixel')
ylabel('Bias')
subplot(3,1,3)
plot(dc)
xlabel('Pixel')
ylabel('DC (counts/s)')

% saves bias and dark_current as files.    
save(biasfile,'bias','-ASCII')
save(dcfile,'dc','-ASCII')
    
% finds the average temperature.    
    disp(['Mean CCD TEMP: ' num2str(mean(Tmean))]);
    disp(['Min CCD TEMP: ' num2str(min(Tmean))]);
    disp(['Max CCD TEMP: ' num2str(max(Tmean))]);