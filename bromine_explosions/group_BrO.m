function [ out_arr, out_arr2 ] = group_BrO(bee_dataset, run_start, run_end, type, varargin)
% [ out_arr ] = group_BrO( run_start, run_end, type, params, threshold)
% 
% Function to find dominant parameters in chunks of the BrO dataset
% INPUT:
%   run_start, run_end: datetime arrays specifying the start and end of the
%                       time periods of interest
%   type: string of output that's needed, extra inputs detailed below:
%
%       'wdir': gives dominant wind direction index in each time period. 
%               Next argument must be the threshold to detemine whech index
%               is dominant. Threshold must be greater than 0.5.
%               Returns:
%                   1: northerly winds dominate
%                   2: southeasterly winds dominate
%                   3: other wind directions dominate
%                   0: no index exceeds the dominance threshold
%                   4: number of missing datapoints exceeds the dominance threshold
%                   Also returns mean wind dir calculated from raw PWS data,
%                   were indices 1, 2, 3 are the same, and 4 indicates that
%                   wdir data is missing for the entire period
%                   (different since wdir for each BrO entry is already
%                   grouped, so deault output groups wdir twice, while
%                   second output groups it only once)
%
%       'bro': groups BrO columns by the bins given in the next variable.
%              The input should be the upper bounds of the bins, and the
%              last bin is considered to go to infinity.
%              Output is the index of the bin the mean BrO column falls
%              into for the given period.
%       'ssa': same as 'bro', for OPC supermicron particle concentrations
%       'o3': same as 'bro', for surface ozone vmr (in ppb)
%       'inv': same as 'bro', for surface to lab inversion strength
%
% OUTPUT:
%   out_arr: array of indices (size of run start and run end), as specified
%            under input 'type' 
%   out_arr2: mean values for requested output (slightly different for 'wdir') 
%
% @Kristof Bognar, 2019


% output variable
out_arr=NaN(length(run_start),1);

% set later where relevant
do_sort=0;

% if alternate output is requested
if nargout==2
    out_arr2=NaN(length(run_start),1);
    
    % if alternate output is mean wdir
    if strcmp(type,'wdir')
        load('/home/kristof/work/weather_stations/ridge_lab/PWS_all.mat');
        data((month(data.DateTime)>5 | month(data.DateTime)<3),:)=[];
        pws_data=data;
    end
end

for i=1:length(run_start)
    
    % get data within current time interval
    % use equal on both ends since some times are at either end of
    % measurement gaps, and both need to be included
    tmp=bee_dataset(bee_dataset.times>=run_start(i) & bee_dataset.times<=run_end(i),:);
    
    switch type
        case 'wdir' 
            %%
            
            % check input
            if varargin{1}<0.5 || varargin{1}>1
                error('Need 0.5 < threshold <= 1');
            end
            
            % replace NaNs
            tmp.N_SE_rest(isnan(tmp.N_SE_rest))=4;
            
            % find how many times each index repeats
            [repeats,values]=hist(tmp.N_SE_rest,unique(tmp.N_SE_rest));     

            try
                % if one index exceeds dominance threshold, that's the output
                out_arr(i)=values(repeats>=size(tmp,1)*varargin{1});
            catch
                % if none are dominant, output is 0
                out_arr(i)=0;
            end
            
            % alternative approach: use 'find_coincident_mean.m' with time
            % periods of interest
            if nargout==2
                out_arr2(i)=find_coincident_mean(mean([run_start(i),run_end(i)]),...
                                                 pws_data.DateTime, pws_data.WindDir,...
                                                 minutes(run_end(i)-run_start(i))/2, true);
                if isnan(out_arr2(i))
                    out_arr2(i)=4;
                elseif out_arr2(i) >= 324 || out_arr2(i) < 24
                    out_arr2(i)=1;
                elseif out_arr2(i) >= 93 && out_arr2(i) < 153
                    out_arr2(i)=2;
                else
                    out_arr2(i)=3;
                end
            end
            
        case 'bro'
            % get mean column in selected period
            mean_tmp=mean(tmp.bro_col);
            do_sort=1;
            
        case 'ssa_hm'
            % get mean aer concetration for sizes half micron and up
            mean_tmp=mean(tmp.aer_halfmicron);
            do_sort=1;

        case 'ssa_sm'
            % get mean aer concetration for supermicron particles
            mean_tmp=mean(tmp.aer_supermicron);
            do_sort=1;
            
        case 'o3'
            % get mean surface ozone
            mean_tmp=mean(tmp.o3_surf);
            do_sort=1;

        case 'inv'
            % get mean T inversion between surface and lab
            mean_tmp=mean(tmp.sonde_dT);
            do_sort=1;
            
    end
    
    if do_sort
        
        % sort selected mean into the provided list of upper bounds
        [~,b]=sort([varargin{1}, mean_tmp]);

        % location of largest index is the bin the mean value falls
        % into, equivalent to mean >= bin edge
        [~,out_arr(i)]=max(b);
        
        % save mean value as well
        if nargout==2, out_arr2(i)=mean_tmp; end
            
    end
    
end

end

