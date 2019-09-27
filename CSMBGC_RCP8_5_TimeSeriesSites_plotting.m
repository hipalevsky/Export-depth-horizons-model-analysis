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
    end
end

%%%%% Next step is to calculate e-ratios and visualize output thus far

%%%% Also do some testing of Taylor decomposition method intended to employ
%%%% globally just at these four time series stations