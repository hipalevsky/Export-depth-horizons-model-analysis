%%% Load output extracted in CSMBGC_RCP8_5_TimeSeriesSites.m
load RCP8_5_2005to2100_TimeSeriesSites.mat
    nstn = length(latstn);
    len = length(z);
    [len_NPP,~,~] = size(out_ts{1}.NPP_diat);
    nummon = 12*length(yrslist);

%% Extract maximum annual MLD
for i = 1:nstn
    out_ts{i}.mld_max = max(out_ts{i}.mld);
end

%% Determine maximum annual MLD from beginning of the century
% Do this for both 10 and 20 yr windows to check sensitivity
% Goal is to use this as fixed depth criteria for full century to enable
% Taylor decomposition
for i = 1:nstn
    out_ts{i}.mld_max_beg10yrs = mean(out_ts{i}.mld_max(1:10));
    out_ts{i}.mld_max_beg20yrs = mean(out_ts{i}.mld_max(1:20));
end

%% Map maximum annual MLD onto out_ts_plot shape
for i = 1:length(yrslist)
    for j = 1:nstn
        out_ts_plot{j}.mldmax((i-1)*12+1:i*12) = out_ts{j}.mld_max(i);
    end
end

%% Create extended NPP arrays
for j = 1:nstn
    out_ts_plot{j}.NPP_diat_ext = repmat(out_ts_plot{j}.NPP_diat(end,:), len, 1);
        out_ts_plot{j}.NPP_diat_ext(1:len_NPP,:) = out_ts_plot{j}.NPP_diat;
    out_ts_plot{j}.NPP_diaz_ext = repmat(out_ts_plot{j}.NPP_diaz(end,:), len, 1);
        out_ts_plot{j}.NPP_diaz_ext(1:len_NPP,:) = out_ts_plot{j}.NPP_diaz;
    out_ts_plot{j}.NPP_sp_ext = repmat(out_ts_plot{j}.NPP_sp(end,:), len, 1);
        out_ts_plot{j}.NPP_sp_ext(1:len_NPP,:) = out_ts_plot{j}.NPP_sp;
    out_ts_plot{j}.NPP_sum_ext = repmat(out_ts_plot{j}.NPP_sum(end,:), len, 1);
        out_ts_plot{j}.NPP_sum_ext(1:len_NPP,:) = out_ts_plot{j}.NPP_sum;
end

%% Determine POCflux and NPP values at mldmax
for i = 1:nummon
    for j = 1:nstn
        %Calculate values at mldmax
        out_ts_plot{j}.POCflux_mldmax(i) = interp1(z_top, squeeze(out_ts_plot{j}.POCflux(:,i)), out_ts_plot{j}.mldmax(i));
        out_ts_plot{j}.NPP_diat_mldmax(i) = interp1(z_top, squeeze(out_ts_plot{j}.NPP_diat_ext(:,i)), out_ts_plot{j}.mldmax(i));
        out_ts_plot{j}.NPP_diaz_mldmax(i) = interp1(z_top, squeeze(out_ts_plot{j}.NPP_diaz_ext(:,i)), out_ts_plot{j}.mldmax(i));
        out_ts_plot{j}.NPP_sp_mldmax(i) = interp1(z_top, squeeze(out_ts_plot{j}.NPP_sp_ext(:,i)), out_ts_plot{j}.mldmax(i));
        out_ts_plot{j}.NPP_sum_mldmax(i) = interp1(z_top, squeeze(out_ts_plot{j}.NPP_sum_ext(:,i)), out_ts_plot{j}.mldmax(i));
        %Calculate values at beginning of century mldmax (just POCflux and NPP_sum)
        out_ts_plot{j}.POCflux_mldmax_beg10yrs(i) = interp1(z_top, squeeze(out_ts_plot{j}.POCflux(:,i)), out_ts{j}.mld_max_beg10yrs);
        out_ts_plot{j}.POCflux_mldmax_beg20yrs(i) = interp1(z_top, squeeze(out_ts_plot{j}.POCflux(:,i)), out_ts{j}.mld_max_beg20yrs);
        out_ts_plot{j}.NPP_sum_mldmax_beg10yrs(i) = interp1(z_top, squeeze(out_ts_plot{j}.NPP_sum_ext(:,i)), out_ts{j}.mld_max_beg10yrs);
        out_ts_plot{j}.NPP_sum_mldmax_beg20yrs(i) = interp1(z_top, squeeze(out_ts_plot{j}.NPP_sum_ext(:,i)), out_ts{j}.mld_max_beg20yrs);
    end
end
%% Calculate annual mean POC flux and NPP values at mldmax and at 100m
for i = 1:length(yrslist)
    for j = 1:nstn
      %Calculate annual mean values
        out_ts_plot{j}.POCflux_mldmax_annmean(i) = mean(out_ts_plot{j}.POCflux_mldmax((i-1)*12+1:i*12));
        out_ts_plot{j}.POCflux_100m_annmean(i) = mean(out_ts_plot{j}.POCflux(11,(i-1)*12+1:i*12));
        out_ts_plot{j}.NPP_diat_mldmax_annmean(i) = mean(out_ts_plot{j}.NPP_diat_mldmax((i-1)*12+1:i*12));
        out_ts_plot{j}.NPP_diat_100m_annmean(i) = mean(out_ts_plot{j}.NPP_diat(11,(i-1)*12+1:i*12));
        out_ts_plot{j}.NPP_diaz_mldmax_annmean(i) = mean(out_ts_plot{j}.NPP_diaz_mldmax((i-1)*12+1:i*12));
        out_ts_plot{j}.NPP_diaz_100m_annmean(i) = mean(out_ts_plot{j}.NPP_diaz(11,(i-1)*12+1:i*12));
        out_ts_plot{j}.NPP_sp_mldmax_annmean(i) = mean(out_ts_plot{j}.NPP_sp_mldmax((i-1)*12+1:i*12));
        out_ts_plot{j}.NPP_sp_100m_annmean(i) = mean(out_ts_plot{j}.NPP_sp(11,(i-1)*12+1:i*12));
        out_ts_plot{j}.NPP_sum_mldmax_annmean(i) = mean(out_ts_plot{j}.NPP_sum_mldmax((i-1)*12+1:i*12));
        out_ts_plot{j}.NPP_sum_100m_annmean(i) = mean(out_ts_plot{j}.NPP_sum(11,(i-1)*12+1:i*12));
      %Values for mldmax from beginning of century
        out_ts_plot{j}.POCflux_mldmax_beg10yrs_annmean(i) = mean(out_ts_plot{j}.POCflux_mldmax_beg10yrs((i-1)*12+1:i*12));
        out_ts_plot{j}.POCflux_mldmax_beg20yrs_annmean(i) = mean(out_ts_plot{j}.POCflux_mldmax_beg20yrs((i-1)*12+1:i*12));
        out_ts_plot{j}.NPP_sum_mldmax_beg10yrs_annmean(i) = mean(out_ts_plot{j}.NPP_sum_mldmax_beg10yrs((i-1)*12+1:i*12));
        out_ts_plot{j}.NPP_sum_mldmax_beg20yrs_annmean(i) = mean(out_ts_plot{j}.NPP_sum_mldmax_beg20yrs((i-1)*12+1:i*12));
      %Profile values as annual mean
        out_ts_plot{j}.POCflux_annmean(:,i) = mean(out_ts_plot{j}.POCflux(:,(i-1)*12+1:i*12)');
        out_ts_plot{j}.NPP_sum_annmean(:,i) = mean(out_ts_plot{j}.NPP_sum_ext(:,(i-1)*12+1:i*12)');
    end
end

%% Visualize output for each time series station
timeinmonths = [yrslist(1)-1:(1/12):yrslist(end)-1/12]; %discretize to months for plotting
yrstoavg = 10; %years to use in running average

for j = 1:4
    figure(j); clf
        subplot(4,1,1)
    plot(timeinmonths, out_ts_plot{j}.mld,'b-'); hold on;
    plot(yrslist, out_ts_plot{j}.mldmax(1:12:end),'k.','markersize',15); hold on;
    %errorbar(yrslist, movmean(out_ts_plot{j}.mldmax(1:12:end), yrstoavg), movstd(out_ts_plot{j}.mldmax(1:12:end), yrstoavg),'k--'); hold on;
    plot(yrslist, movmean(out_ts_plot{j}.mldmax(1:12:end), yrstoavg),'k--'); hold on;
    axis ij
    legend('MLD','MLD_{annual max}','10-yr mean, MLD_{annual max}','location','northwest')
    title(stnname{j})
        subplot(4,1,2)
    plot(yrslist, out_ts_plot{j}.POCflux_mldmax_annmean,'k-'); hold on;
    plot(yrslist, out_ts_plot{j}.POCflux_100m_annmean,'b-'); hold on;
    legend('Annual POC flux_{MLDmax}','Annual POC flux_{100 m}','location','northwest')
        subplot(4,1,3)
    %plot(yrslist, out_ts_plot{j}.NPP_sum_mldmax_annmean,'k-'); hold on;
    plot(yrslist, out_ts_plot{j}.NPP_sum_100m_annmean,'b-'); hold on;
    plot(yrslist, out_ts_plot{j}.NPP_sp_100m_annmean,'m-'); hold on;
    plot(yrslist, out_ts_plot{j}.NPP_diat_100m_annmean,'g-'); hold on;
    %plot(yrslist, out_ts_plot{j}.NPP_diaz_100m_annmean,'y-'); hold on;
    legend('Annual NPP_{100 m}','Annual NPP_{100 m} from small phytos','Annual NPP_{100 m} from diatoms','location','northwest')
        subplot(4,1,4)
    plot(yrslist, out_ts_plot{j}.POCflux_mldmax_annmean./out_ts_plot{j}.NPP_sum_mldmax_annmean,'k-'); hold on;
    plot(yrslist, out_ts_plot{j}.POCflux_100m_annmean./out_ts_plot{j}.NPP_sum_100m_annmean,'b-'); hold on;
    legend('Annual e-ratio_{MLDmax}','Annual e-ratio_{100 m}','location','northwest')
end

%% Make depth vs. time contour plots of POC flux and e-ratio across all time-series sites
figure(5); clf
M = 10; L = 1.5;
for j = 1:4
        subplot(4,2,j*2-1)
    contourf(yrslist, z_top, out_ts_plot{j}.POCflux_annmean,'linecolor','none'); hold on;
    colormap(cmocean('algae')); h = colorbar; %caxis([0 5]); 
    ylim([0 max(out_ts_plot{j}.mldmax)*1.1])
    plot(yrslist, out_ts_plot{j}.mldmax(1:12:end),'k.','markersize',M); hold on;
    plot(yrslist, movmean(out_ts_plot{j}.mldmax(1:12:end), yrstoavg),'k-','linewidth',L); hold on;
    axis ij
    title([stnname{j} ' POC flux'])
    ylabel('Depth (m)')
    ylabel(h, 'mol C m^{-2} yr^{-1}')
                subplot(4,2,j*2)
    contourf(yrslist, z_top, out_ts_plot{j}.POCflux_annmean./out_ts_plot{j}.NPP_sum_annmean,'linecolor','none'); hold on;
    colormap(cmocean('algae')); colorbar; %caxis([0 0.35]);
    ylim([0 max(out_ts_plot{j}.mldmax)*1.1])
    plot(yrslist, out_ts_plot{j}.mldmax(1:12:end),'k.','markersize',M); hold on;
    plot(yrslist, movmean(out_ts_plot{j}.mldmax(1:12:end), yrstoavg),'k-','linewidth',L); hold on;
    axis ij
    title([stnname{j} ' e-ratio'])
    ylabel('Depth (m)')
%             subplot(4,2,j*2)
%     contourf(yrslist, z_top, out_ts_plot{j}.NPP_sum_annmean,'linecolor','none'); hold on;
%     colormap(cmocean('algae')); colorbar; caxis([0 15]);
%     ylim([0 max(out_ts_plot{j}.mldmax)*1.1])
%     plot(yrslist, out_ts_plot{j}.mldmax(1:12:end),'k.','markersize',15); hold on;
%     plot(yrslist, movmean(out_ts_plot{j}.mldmax(1:12:end), yrstoavg),'k--'); hold on;
%     axis ij
%     title([stnname{j} ' NPP'])
end

%% Calculate Taylor decompositon where:
% dEP/dt = (dNPP/dt x e-ratio)_beginningcenturyMLDmax + (dEratio/dt x NPP)_beginningcenturyMLDmax
%    + (dPOCflux/dz x dMLDmaxdt) + Residual

for j = 1:4
%Calculate Taylor decomposition terms for each time series site
    dEPdt_10yr(j) = mean(out_ts_plot{j}.POCflux_mldmax_annmean(end-9:end)) - mean(out_ts_plot{j}.POCflux_mldmax_annmean(1:10));
    dEPdt_20yr(j) = mean(out_ts_plot{j}.POCflux_mldmax_annmean(end-19:end)) - mean(out_ts_plot{j}.POCflux_mldmax_annmean(1:20));

    dNPPdt_10yr(j) = mean(out_ts_plot{j}.NPP_sum_mldmax_beg10yrs_annmean(end-9:end)) - mean(out_ts_plot{j}.NPP_sum_mldmax_beg10yrs_annmean(1:10));
    dNPPdt_20yr(j) = mean(out_ts_plot{j}.NPP_sum_mldmax_beg20yrs_annmean(end-19:end)) - mean(out_ts_plot{j}.NPP_sum_mldmax_beg20yrs_annmean(1:20));

    TaylorNPPterm_10yr(j) = dNPPdt_10yr(j)*...
        mean((out_ts_plot{j}.POCflux_mldmax_beg10yrs_annmean(1:10)./out_ts_plot{j}.NPP_sum_mldmax_beg10yrs_annmean(1:10)));
    TaylorNPPterm_20yr(j) = dNPPdt_20yr(j)*...
        mean((out_ts_plot{j}.POCflux_mldmax_beg20yrs_annmean(1:20)./out_ts_plot{j}.NPP_sum_mldmax_beg20yrs_annmean(1:20)));

    dEratiodt_10yr(j) = mean((out_ts_plot{j}.POCflux_mldmax_beg10yrs_annmean(end-9:end))./(out_ts_plot{j}.NPP_sum_mldmax_beg10yrs_annmean(end-9:end))) -...
        mean((out_ts_plot{j}.POCflux_mldmax_beg10yrs_annmean(1:10))./(out_ts_plot{j}.NPP_sum_mldmax_beg10yrs_annmean(1:10)));
    dEratiodt_20yr(j) = mean((out_ts_plot{j}.POCflux_mldmax_beg10yrs_annmean(end-19:end))./(out_ts_plot{j}.NPP_sum_mldmax_beg10yrs_annmean(end-19:end))) -...
        mean((out_ts_plot{j}.POCflux_mldmax_beg10yrs_annmean(1:20))./(out_ts_plot{j}.NPP_sum_mldmax_beg10yrs_annmean(1:20)));

    TaylorEratioterm_10yr(j) = dEratiodt_10yr(j)*mean(out_ts_plot{j}.NPP_sum_mldmax_beg10yrs_annmean(1:10));
    TaylorEratioterm_20yr(j) = dEratiodt_20yr(j)*mean(out_ts_plot{j}.NPP_sum_mldmax_beg20yrs_annmean(1:20));

    TaylorMLDmaxterm_10yr(j) = mean(out_ts_plot{j}.POCflux_mldmax_annmean(end-9:end)) - mean(out_ts_plot{j}.POCflux_mldmax_beg10yrs_annmean(end-9:end));
    TaylorMLDmaxterm_20yr(j) = mean(out_ts_plot{j}.POCflux_mldmax_annmean(end-19:end)) - mean(out_ts_plot{j}.POCflux_mldmax_beg20yrs_annmean(end-19:end));

    TaylorResidual_10yr(j) = dEPdt_10yr(j) - TaylorNPPterm_10yr(j) - TaylorEratioterm_10yr(j) - TaylorMLDmaxterm_10yr(j);
    TaylorResidual_20yr(j) = dEPdt_20yr(j) - TaylorNPPterm_20yr(j) - TaylorEratioterm_20yr(j) - TaylorMLDmaxterm_20yr(j);
    
    NormalizeData_10yr(j) = mean(out_ts_plot{j}.POCflux_mldmax_annmean(1:10));
    NormalizeData_20yr(j) = mean(out_ts_plot{j}.POCflux_mldmax_annmean(1:20));
    
end

%% Plot bar graph with results of Taylor decomposition
figure(6); clf
for i = 1:4
subplot(4, 2, i*2-1)
    bar([dEPdt_10yr(i) TaylorNPPterm_10yr(i) TaylorEratioterm_10yr(i) TaylorMLDmaxterm_10yr(i) TaylorResidual_10yr(i)])
    title([stnname{i} ' 10-yr Taylor decomposition'])
    xticklabels([{'dEP/dt','dNPP/dt','dE-ratio/dt','dMLDmax/dt','Residual'}])
subplot(4, 2, i*2)
    bar([dEPdt_20yr(i) TaylorNPPterm_20yr(i) TaylorEratioterm_20yr(i) TaylorMLDmaxterm_20yr(i) TaylorResidual_20yr(i)])
    title([stnname{i} ' 20-yr Taylor decomposition'])
    xticklabels([{'dEP/dt','dNPP/dt','dE-ratio/dt','dMLDmax/dt','Residual'}])
end

%% Plot a prettier version of the bar graph with results of Taylor decomposition
figure(7); clf
    bar([dEPdt_20yr; TaylorNPPterm_20yr; TaylorEratioterm_20yr; TaylorMLDmaxterm_20yr; TaylorResidual_20yr])
    title(['Taylor decomposition of POC flux changes from 2005-2021 to 2081-2100 at stations of interest'])
    xticklabels([{'Change in POCflux_{MLDmax}','NPP component','e-ratio component','\DeltaMLD_{max} component','Residual'}])
    legend(stnname)
    ylabel('mol C m^{-2} yr^{-1}')


