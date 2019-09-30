function sp = get_spectrum_v2(filename,start_time)
% June 2003, by Jennifer Walker
% Returns the vector sp(1:2000) which is the spectrum for the cycle from start_time to end_time.
% The arguments start_time and end_time should be strings, written as the time appears in the
% file, e.g. start_time = '3:38:12'.  The directory containing the input file is
% the variable specdir, which is assigned at the beginning of this function.
%
% MODIFIED FEB. 2005 FOR NEW CCD WITH 2048 PIXELS!!  - Annemarie Fraser
% Modified Oct, 2010 for new file-naming process - Cristen Adams

filepath = filename;

fid = fopen(filepath,'r');
if fid == -1 	% file was not opened
    error(strcat('Error opening file [',filepath,'] - file not opened.'));
end

fgetl(fid);     %Skip header lines
fgetl(fid);

while ~feof(fid)
    str=fgetl(fid);
    disp(str);
    if findstr(str,start_time)
        disp('found start time');
        for i=1:2048
            sp(i)=fscanf(fid,'%f',1);
        end
        found =1;
        fclose(fid);
        return;
    end
    for i=1:2048
        fgetl(fid);
    end
end
fclose(fid);
error('End of file reached and start_time/end_time not found');
