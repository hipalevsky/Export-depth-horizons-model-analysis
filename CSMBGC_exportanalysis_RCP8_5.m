% Analysis of export at varying depth horions over 21st century in RCP 8.5
% CESM output

% Initial goal is to pull out climatological data binned over 1st and last
% 20 years of output (2005-2024 and 2081-2100). Start by writing code to
% just do for initial 20 year period.

%Specify size of output to be extracted
d1s = 320; d2s = 384; %lat-lon grid size
zlen = 60; zlen_NPP = 15; %depth grid sizes

%Choose years to include in climatology
yrslist = [2005:2006];

%Initialize arrays to hold output from all years
mld_all = NaN*ones(d1s,d2s,12*length(yrslist));
POCflux_all = NaN*ones(d1s,d2s,zlen,12*length(yrslist));
POCprod_all = NaN*ones(d1s,d2s,zlen,12*length(yrslist));
NPP_diat_all = NaN*ones(d1s,d2s,zlen_NPP,12*length(yrslist));
NPP_diaz_all = NaN*ones(d1s,d2s,zlen_NPP,12*length(yrslist));
NPP_sp_all = NaN*ones(d1s,d2s,zlen_NPP,12*length(yrslist));
NPP_sum_all = NaN*ones(d1s,d2s,zlen_NPP,12*length(yrslist));

%Extract all data from years to use in making climatology
for i = 1:length(yrslist)
    [z, z_top, TLAT, TLONG, TAREA, REGION_MASK, mld_yr, POCflux_yr, POCprod_yr, NPP_diat_yr, NPP_diaz_yr, NPP_sp_yr, NPP_sum_yr] = CSMBGC_openfiles_RCP8_5(yrslist(i)); 
    mld_all(:,:,(1+12*(i-1)):(12+12*(i-1))) = mld_yr;
    POCflux_all(:,:,:,(1+12*(i-1)):(12+12*(i-1))) = POCflux_yr;
    POCprod_all(:,:,:,(1+12*(i-1)):(12+12*(i-1))) = POCprod_yr;
    NPP_diat_all(:,:,:,(1+12*(i-1)):(12+12*(i-1))) = NPP_diat_yr;
    NPP_diaz_all(:,:,:,(1+12*(i-1)):(12+12*(i-1))) = NPP_diaz_yr;
    NPP_sp_all(:,:,:,(1+12*(i-1)):(12+12*(i-1))) = NPP_sp_yr;
    NPP_sum_all(:,:,:,(1+12*(i-1)):(12+12*(i-1))) = NPP_sum_yr;
    disp({'Year ' num2str(yrslist(i)) ' extracted at ' datestr(now,13)})
end


%% Calculate depth criteria - TRY WITHOUT INTERPOLATING TO FINER DEPTH GRID BECAUSE THAT TAKES A LONG TIME TO RUN
%Euphotic depth, as defined by Lima et al. 2014: depth where POC production is equal to 1% of maximum POC production in the water column
[zeu] = zeuFromPOCProd(POCprod_all, z);

%Particle compensation depth: depth of maximum POC flux
[zcomp, POCflux_zcomp] = zcompFromPOCFlux(z_top, POCflux_all);

%Calculate maximum annual MLD from seasonal MLD
mldmaxyr = NaN*ones(d1s,d2s,length(yrslist)); %initialize array with to hold mldmax for each year
mldmax = NaN*ones(d1s,d2s,12*length(yrslist)); %initialize array with time same size as input
for i = 1:length(yrslist)
    mldmaxyr(:,:,i) = max(mld_all(:,:,(i-1)*12+1:i*12),[],3);
    mldmax(:,:,(i-1)*12+1:i*12) = repmat(mldmaxyr(:,:,i),1,1,12);
end

disp({'Depth criteria calculated at ' datestr(now,13)})

%% Calculate values at depth criteria depths - NOTE THAT I SWITCHED ORDER
%Create a cell array with all the depth criteria. 
    depthCriteria{1} = mld_all;
    depthCriteria{2} = zcomp; 
    depthCriteria{3} = zeu; 
    depthCriteria{4} = mldmax;
    depthCriteria{5} = 100*ones(d1s,d2s,12*length(yrslist));
POCflux_depthCriteria = depthCriteriaVals(POCflux_all, z_top, depthCriteria);

disp({'POC flux at depth criteria calculated at ' datestr(now,13)})

%% Create NPP extended that assumes no additional NPP below 145 m (but doesn't fill with NaNs)
    NPP_sum_extended = NaN*ones(d1s,d2s,length(z),12*length(yrslist));
        NPP_sum_extended(:,:,1:15,:) = NPP_sum_all;
        NPP_sum_extended(:,:,16:end,:) = repmat(NPP_sum_all(:,:,15,:),[1,1,45,1]);
    NPP_diat_extended = NaN*ones(d1s,d2s,length(z),12*length(yrslist));
        NPP_diat_extended(:,:,1:15,:) = NPP_diat_all;
        NPP_diat_extended(:,:,16:end,:) = repmat(NPP_diat_all(:,:,15,:),[1,1,45,1]);
    NPP_diaz_extended = NaN*ones(d1s,d2s,length(z),12*length(yrslist));
        NPP_diaz_extended(:,:,1:15,:) = NPP_diaz_all;
        NPP_diaz_extended(:,:,16:end,:) = repmat(NPP_diaz_all(:,:,15,:),[1,1,45,1]);
    NPP_sp_extended = NaN*ones(d1s,d2s,length(z),12*length(yrslist));
        NPP_sp_extended(:,:,1:15,:) = NPP_sp_all;
        NPP_sp_extended(:,:,16:end,:) = repmat(NPP_sp_all(:,:,15,:),[1,1,45,1]);
    
NPPsum_depthCriteria = depthCriteriaVals(NPP_sum_extended, z, depthCriteria);
    NPPdiat_depthCriteria = depthCriteriaVals(NPP_diat_extended, z, depthCriteria);
    NPPdiaz_depthCriteria = depthCriteriaVals(NPP_diaz_extended, z, depthCriteria);
    NPPsp_depthCriteria = depthCriteriaVals(NPP_sp_extended, z, depthCriteria);

disp({'NPP at depth criteria calculated at ' datestr(now,13)})

%% Calculate climatological (over all years in yrslist) and annual mean values

numCriteria = length(depthCriteria);
year = repmat(yrslist,12,1); year = year(:);

for i = 1:numCriteria
    %%%% For depth criteria
    depthCriteria_clim{i} = climatologicalMean(depthCriteria{i}); % Climatological monthly values (mean over all years)
    depthCriteria_wAnnMean{i} = weightedAnnualMean(depthCriteria{i},... % Annual mean values for each year, weighted by monthly POC flux
        POCflux_depthCriteria{i}, year); 
    depthCriteria_wAnnMean_clim{i} = mean(depthCriteria_wAnnMean{i},3); % Climatological annual POC flux-weighted mean for each grid cell
    %%%% For POC flux rates
    POCflux_depthCriteria_clim{i} = climatologicalMean(POCflux_depthCriteria{i}); % Climatological monthly values (mean over all years)
    for k = 1:length(yrslist) % Annual mean values for each year
        idyr = find(year == yrslist(k));
        POCflux_depthCriteria_AnnMean{i}(:,:,k) = mean(POCflux_depthCriteria{i}(:,:,idyr),3);
    end
    POCflux_depthCriteria_AnnMean_clim{i} = mean(POCflux_depthCriteria_AnnMean{i},3); % Climatological annual POC flux rates for each grid cell at each criterion
    %%%% For NPP rates
    NPPsum_depthCriteria_clim{i} = climatologicalMean(NPPsum_depthCriteria{i}); % Climatological monthly values (mean over all years)
    NPPdiat_depthCriteria_clim{i} = climatologicalMean(NPPdiat_depthCriteria{i}); % Climatological monthly values (mean over all years)
    NPPdiaz_depthCriteria_clim{i} = climatologicalMean(NPPdiaz_depthCriteria{i}); % Climatological monthly values (mean over all years)
    NPPsp_depthCriteria_clim{i} = climatologicalMean(NPPsp_depthCriteria{i}); % Climatological monthly values (mean over all years)
    for k = 1:length(yrslist) % Annual mean values for each year
        idyr = find(year == yrslist(k));
        NPPsum_depthCriteria_AnnMean{i}(:,:,k) = mean(NPPsum_depthCriteria{i}(:,:,idyr),3);
        NPPdiat_depthCriteria_AnnMean{i}(:,:,k) = mean(NPPdiat_depthCriteria{i}(:,:,idyr),3);
        NPPdiaz_depthCriteria_AnnMean{i}(:,:,k) = mean(NPPdiaz_depthCriteria{i}(:,:,idyr),3);
        NPPsp_depthCriteria_AnnMean{i}(:,:,k) = mean(NPPsp_depthCriteria{i}(:,:,idyr),3);
    end
    NPPsum_depthCriteria_AnnMean_clim{i} = mean(NPPsum_depthCriteria_AnnMean{i},3); % Climatological annual NPP rates for each grid cell at each criterion
    NPPdiat_depthCriteria_AnnMean_clim{i} = mean(NPPdiat_depthCriteria_AnnMean{i},3); % Climatological annual NPP rates for each grid cell at each criterion
    NPPdiaz_depthCriteria_AnnMean_clim{i} = mean(NPPdiaz_depthCriteria_AnnMean{i},3); % Climatological annual NPP rates for each grid cell at each criterion
    NPPsp_depthCriteria_AnnMean_clim{i} = mean(NPPsp_depthCriteria_AnnMean{i},3); % Climatological annual NPP rates for each grid cell at each criterion
end

%%%% For POC flux at all depths
POCflux_all_clim = climatologicalMean(POCflux_all);
for k = 1:length(yrslist)
    idyr = find(year == yrslist(k));
    POCflux_all_AnnMean(:,:,:,k) = mean(POCflux_all(:,:,:,idyr),4);
end
POCflux_all_AnnMean_clim = mean(POCflux_all_AnnMean,4);

%%%% For integrated NPP at all depths
NPPsum_all_clim = climatologicalMean(NPP_sum_extended);
for k = 1:length(yrslist)
    idyr = find(year == yrslist(k));
    NPPsum_all_AnnMean(:,:,:,k) = mean(NPP_sum_extended(:,:,:,idyr),4);
end
NPPsum_all_AnnMean_clim = mean(NPPsum_all_AnnMean,4);

disp({'Climatology calculated at ' datestr(now,13)})

%%
save RCP8_5_2005to2006_ClimOnly yrslist depthCriteria_clim depthCriteria_wAnnMean depthCriteria_wAnnMean_clim...
    POCflux_depthCriteria_clim POCflux_depthCriteria_AnnMean POCflux_depthCriteria_AnnMean_clim...
    NPPsum_depthCriteria_clim NPPdiat_depthCriteria_clim NPPdiaz_depthCriteria_clim NPPsp_depthCriteria_clim...
    NPPsum_depthCriteria_AnnMean NPPdiat_depthCriteria_AnnMean NPPdiaz_depthCriteria_AnnMean NPPsp_depthCriteria_AnnMean...
    NPPsum_depthCriteria_AnnMean_clim NPPdiat_depthCriteria_AnnMean_clim NPPdiaz_depthCriteria_AnnMean_clim NPPsp_depthCriteria_AnnMean_clim

save RCP8_5_2005to2006_ClimPOCNPPsections yrslist POCflux_all_clim POCflux_all_AnnMean POCflux_all_AnnMean_clim...
    NPPsum_all_clim NPPsum_all_AnnMean NPPsum_all_AnnMean_clim

%% Visualizing modern data (after the analysis of the Lima 2014 data)
%% Calculate regional and global rates
    latint = 4;
latdiv = [-66:latint:66]; %define latitude divisions for zonal regions

clear NPP_sum_glob NPP_sum_reg POCflux_sum_glob POCflux_sum_reg area_reg
for i = 1:numCriteria
    [NPP_sum_glob(i), NPP_sum_reg(:,i), area_reg] = areaRateSum(NPPsum_depthCriteria_AnnMean_clim{i},latdiv,TLAT,TAREA,REGION_MASK);
    [POCflux_sum_glob(i), POCflux_sum_reg(:,i), area_reg] = areaRateSum(POCflux_depthCriteria_AnnMean_clim{i},latdiv,TLAT,TAREA,REGION_MASK);
end

%% Need to calculate not as rate for depths
for j = 1:numCriteria
    A = depthCriteria_wAnnMean_clim{j};
    for i = 1:length(area_reg)
        reg_id = find(TLAT > latdiv(i) & TLAT <= latdiv(i+1) & REGION_MASK > 0);
        depthCriteria_wAnnMean_clim_regmean(i,j) = nanmean(A(reg_id));
        depthCriteria_wAnnMean_clim_regstd(i,j) = nanstd(A(reg_id));
    end
end

%% Regrid for plotting
%Note that if the regridded grid size gets too fine for the input data,
%this function creates anomalous patterns where input data is scarce
lonmin = 0; lonmax = 360;
latmin = -85; latmax = 85;
lon_gridsize = 2; lat_gridsize = 2;

[glon, glat, POCflux_interp_AnnMean_clim_grid] = regrid_even(TLONG, TLAT, POCflux_all_AnnMean_clim,...
    lonmin, lonmax, latmin, latmax, lon_gridsize, lat_gridsize);
[glon, glat, NPPsum_interp_AnnMean_clim_grid] = regrid_even(TLONG, TLAT, NPPsum_all_AnnMean_clim,...
    lonmin, lonmax, latmin, latmax, lon_gridsize, lat_gridsize);

% Pack depth criteria, POC flux at depth criteria, and NPP sum at depth criteria into stack
stackedCriteria(:,:,1) = min(depthCriteria_clim{1},[],3); % Add minimum annual mixed layer depth at top of stack
for i = 1:numCriteria
    stackedCriteria(:,:,i+1) = depthCriteria_wAnnMean_clim{i};
    stackedCriteria(:,:,i+numCriteria+1) = POCflux_depthCriteria_AnnMean_clim{i};
    stackedCriteria(:,:,i+2*numCriteria+1) = NPPsum_depthCriteria_AnnMean_clim{i};
end

[glon, glat, stackedCriteria_grid] = regrid_even(TLONG, TLAT, stackedCriteria,...
    lonmin, lonmax, latmin, latmax, lon_gridsize, lat_gridsize);

%% Define regions based on relationships among weighted annual mean depth criteria - use depthCriteria_wAnnMean_clim

% Calculate relationships with maximum annual MLD to create regions
[m,n] = size(depthCriteria_wAnnMean_clim{1});
depthregID = zeros(m,n);
for j = 1:m
    for k = 1:n
        %Assign regions based on relationship between maximum annual
        %MLD and other depth criteria
        if depthCriteria_wAnnMean_clim{4}(j,k) > depthCriteria_wAnnMean_clim{2}(j,k)
            depthregID(j,k) = 1; %maximum mld > euphotic depth
        elseif depthCriteria_wAnnMean_clim{4}(j,k) <= depthCriteria_wAnnMean_clim{2}(j,k) & depthCriteria_wAnnMean_clim{4}(j,k) > depthCriteria_wAnnMean_clim{3}(j,k)
            depthregID(j,k) = 2; %euphotic depth > maximum mld > compensation depth
        elseif depthCriteria_wAnnMean_clim{4}(j,k) < depthCriteria_wAnnMean_clim{3}(j,k)
            depthregID(j,k) = 3;
        end
    end
end

%Regrid for plotting
[glon, glat, depthregID_grid] = regrid_even(TLONG, TLAT, depthregID,...
    lonmin, lonmax, latmin, latmax, lon_gridsize, lat_gridsize);

%Time series station locations for plotting
bats = [32 360-64];
hot = [22.75 360-158];
osp = [50 360-145];
keo = [32.3 144.5];
dp1 = [-57 360-64];
dp4 = [-61.5 360-62];
icel = [68 360-12.66];
irming = [59 360-39];
kuroshio = [35 150]; %region 1
eqpac = [0 200]; %region 1
wsubtroppac = [8 360-105]; %region 3
subtropsatl = [-30 360-15]; %region 2.5
southoc_indian = [-40 90];

%These are just the subset of stations currently being used (to plot on maps)
stn_plot = [irming; kuroshio; eqpac; southoc_indian];

stns = [dp1; irming; hot; osp; kuroshio; eqpac; wsubtroppac; subtropsatl; southoc_indian];
stnNames = {'Drake Passage (north)','Irminger Sea','HOT','OSP','Kuroshio','Equatorial Pacific','W. Subtropical N. Pacific','Subtropical S. Atlantic','S. Ocean, Indian'};
zmaxstn = [350 450 150 150 350 150 150 150 450];

%% Run plotting script to make figures
CSMBGC_plotting_RCP8_5
