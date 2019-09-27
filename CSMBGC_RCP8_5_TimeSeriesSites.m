%% Open a single file to extract grid
[z, z_top, TLAT, TLONG, TAREA] = CSMBGC_openfiles_RCP8_5(2005); 

%% Select sites for analysis and identify model grid cells of interest
latstn = [60 57 40 -48];
lonstn = [303 336 170 160];
stnname = {'Labrador Sea','Iceland Basin','Kuroshio Extension','South Pacific'};
    nstn = length(latstn);
for k = 1:nstn
    [inds(k,:),gridpt(k,:)] = findNearestPoint(TLAT,TLONG,latstn(k),lonstn(k));
end

%% Extract output over entire 21st century RCP 8.5 simulation for selected time series sites
yrslist = [2005:2100];
for i = 1:length(yrslist)
    %Extract all output for a given year
    [out] = CSMBGC_openfiles_RCP8_5_AllVars(yrslist(i)); 
    %Pull values just from locations of time series sites
    for j = 1:nstn
        out_ts{j}.TEMP(:,:,i) = squeeze(out.TEMP(inds(j,1),inds(j,2),:,:));
        out_ts{j}.SALT(:,:,i) = squeeze(out.SALT(inds(j,1),inds(j,2),:,:));
        out_ts{j}.ALK(:,:,i) = squeeze(out.ALK(inds(j,1),inds(j,2),:,:));
        out_ts{j}.DIC(:,:,i) = squeeze(out.DIC(inds(j,1),inds(j,2),:,:));
        out_ts{j}.O2(:,:,i) = squeeze(out.O2(inds(j,1),inds(j,2),:,:));
        out_ts{j}.POCflux(:,:,i) = squeeze(out.POCflux(inds(j,1),inds(j,2),:,:));
        out_ts{j}.POCprod(:,:,i) = squeeze(out.POCprod(inds(j,1),inds(j,2),:,:));
        out_ts{j}.NPP_diat(:,:,i) = squeeze(out.NPP_diat(inds(j,1),inds(j,2),:,:));
        out_ts{j}.NPP_diaz(:,:,i) = squeeze(out.NPP_diaz(inds(j,1),inds(j,2),:,:));
        out_ts{j}.NPP_sp(:,:,i) = squeeze(out.NPP_sp(inds(j,1),inds(j,2),:,:));
        out_ts{j}.NPP_sum(:,:,i) = squeeze(out.NPP_sum(inds(j,1),inds(j,2),:,:));
        out_ts{j}.pCO2SURF(:,i) = squeeze(out.pCO2SURF(inds(j,1),inds(j,2),:));
        out_ts{j}.ATM_CO2(:,i) = squeeze(out.ATM_CO2(inds(j,1),inds(j,2),:));
        out_ts{j}.mld(:,i) = squeeze(out.mld(inds(j,1),inds(j,2),:));
        out_ts{j}.F_O2(:,i) = squeeze(out.F_O2(inds(j,1),inds(j,2),:));
        out_ts{j}.F_CO2(:,i) = squeeze(out.F_CO2(inds(j,1),inds(j,2),:));
    end
    clear out
end

%% Reshape data for plotting
len = length(z);
[len_NPP,~,~] = size(out_ts{1}.NPP_diat);
nummon = 12*length(yrslist);
for j = 1:nstn
    out_ts_plot{j}.TEMP = reshape(out_ts{j}.TEMP, len, nummon);
    out_ts_plot{j}.SALT = reshape(out_ts{j}.SALT, len, nummon);
    out_ts_plot{j}.ALK = reshape(out_ts{j}.ALK, len, nummon);
    out_ts_plot{j}.DIC = reshape(out_ts{j}.DIC, len, nummon);
    out_ts_plot{j}.O2 = reshape(out_ts{j}.O2, len, nummon);
    out_ts_plot{j}.POCflux = reshape(out_ts{j}.POCflux, len, nummon);
    out_ts_plot{j}.POCprod = reshape(out_ts{j}.POCprod, len, nummon);
    out_ts_plot{j}.NPP_diat = reshape(out_ts{j}.NPP_diat, len_NPP, nummon);
    out_ts_plot{j}.NPP_diaz = reshape(out_ts{j}.NPP_diaz, len_NPP, nummon);
    out_ts_plot{j}.NPP_sp = reshape(out_ts{j}.NPP_sp, len_NPP, nummon);
    out_ts_plot{j}.NPP_sum = reshape(out_ts{j}.NPP_sum, len_NPP, nummon);
    out_ts_plot{j}.pCO2SURF = reshape(out_ts{j}.pCO2SURF, 1, nummon);
    out_ts_plot{j}.ATM_CO2 = reshape(out_ts{j}.ATM_CO2, 1, nummon);
    out_ts_plot{j}.mld = reshape(out_ts{j}.mld, 1, nummon);
    out_ts_plot{j}.F_O2 = reshape(out_ts{j}.F_O2, 1, nummon);
    out_ts_plot{j}.F_CO2 = reshape(out_ts{j}.F_CO2, 1, nummon);
end

%% Visualize time series data for each variable
L = 2.5
bottomplot = -1500;
for i = 1:nstn
    figure(i); clf
    subplot(211)
        contourf([1:nummon],-z,(out_ts_plot{i}.TEMP),'linecolor','none'); hold on;
        plot([1:nummon],-out_ts_plot{i}.mld,'k-','linewidth',L); hold on;
        ylim([bottomplot 0])
    subplot(212)
        contourf([1:nummon],-z,(out_ts_plot{i}.POCflux),'linecolor','none'); hold on;
        plot([1:nummon],-out_ts_plot{i}.mld,'k-','linewidth',L); hold on;
        ylim([bottomplot 0])
end

%% Save output for further analysis
save RCP8_5_2005to2100_TimeSeriesSites out_ts out_ts_plot latstn lonstn stnname z



