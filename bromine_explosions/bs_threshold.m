function [ output_args ] = bs_threshold()
%BS_THRESHOLD blowing snow threshold estimates


load('/home/kristof/work/BEEs/BEE_dataset_all.mat')
bee_dataset(bee_dataset.times.Year==2015,:)=[];


bs_pws=9.43+0.18*bee_dataset.T_PWS+0.0033*bee_dataset.T_PWS.^2;
bs_ews=9.43+0.18*bee_dataset.T_EWS+0.0033*bee_dataset.T_EWS.^2;


ind_pwd_above=(bee_dataset.wspd_ms>bs_pws & ~isnan(bee_dataset.aer_halfmicron));
ind_pwd_below=(bee_dataset.wspd_ms<=bs_pws & ~isnan(bee_dataset.aer_halfmicron));

disp(['Coarse mode above BS threshold: ' ...
      num2str(round(mean(bee_dataset.aer_halfmicron(ind_pwd_above)),2)) ' +- ' ...
      num2str(round(std(bee_dataset.aer_halfmicron(ind_pwd_above)),2)) ' (' ...
      num2str(round(std(bee_dataset.aer_halfmicron(ind_pwd_above)) / ...
                    sqrt(sum(ind_pwd_above)),2)) ')'])
  
disp(['Coarse mode below BS threshold: ' ...
      num2str(round(mean(bee_dataset.aer_halfmicron(ind_pwd_below)),2)) ' +- ' ...
      num2str(round(std(bee_dataset.aer_halfmicron(ind_pwd_below)),2)) ' (' ...
      num2str(round(std(bee_dataset.aer_halfmicron(ind_pwd_below)) / ...
                    sqrt(sum(ind_pwd_below)),2)) ')'])

% figure()
% histogram(bee_dataset.aer_halfmicron(ind_pwd_below),'normalization','probability'), hold on
% histogram(bee_dataset.aer_halfmicron(ind_pwd_above),'normalization','probability')

% Welch's t-test
[~,p] = ttest2(bee_dataset.aer_halfmicron(ind_pwd_below),...
               bee_dataset.aer_halfmicron(ind_pwd_above),...
               'tail','left','vartype','unequal')

% Mann–Whitney U test
p = ranksum(bee_dataset.aer_halfmicron(ind_pwd_below),...
            bee_dataset.aer_halfmicron(ind_pwd_above),...
            'tail','left')

end

