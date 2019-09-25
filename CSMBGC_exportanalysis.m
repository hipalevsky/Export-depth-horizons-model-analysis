%% Main script to analyze CSMGBC data
% Current script focuses on POC data, but can be built on to incorporate
% NPP, Si flux, CaCO3 flux, etc.

%% Open monthly files from CSMBGC run (last 10 years of Lima et al 2014)
CSMBGC_openfiles

%% Option to visualize the original CSM grid
plotCSMgrid

%% Convert units
    secinyr = 60*60*24*365;
year = floor((time-1)/365); %the "-1" keeps last day of Dec in current yr
    yrslist = unique(year);
    day = (time-1) - year*365 + 1;
z = z_t/100; %depth from surface to midpoint of layer in meters
z_top = z_w/100; %depth from surface to top of layer in meters
mld = HMXL/100; %mixed layer depth in meters

% POC data (flux, production, and remineralization)
POCflux = POC_FLUX_IN.*secinyr/100/1000; %mmol POC/m^3 cm/sec to mol C m-2 yr-1
POCprod = POC_PROD.*secinyr/1000; %mmol POC/m^3/sec to mol C m-3 yr-1
POCremin = POC_REMIN.*secinyr/1000; %mmol POC/m^3/sec to mol C m-3 yr-1
    [d1s,d2s,d3s,d4s] = size(POCflux); %size of dimensions

% NPP data (diatom, diazotroph, and small phytoplankton groups)
    % Volumetric rates
NPPv_diat = photoC_diat.*secinyr/1000; %mmol C/m^3/sec to mol C m-3 yr-1
NPPv_diaz = photoC_diaz.*secinyr/1000; %mmol C/m^3/sec to mol C m-3 yr-1
NPPv_sp = photoC_sp.*secinyr/1000; %mmol C/m^3/sec to mol C m-3 yr-1
    % Depth-integrated rates (note: calculated before interpolation to
    % avoid artefacts with NaNs in surface interpolation)
    dz_grid = shiftdim(repmat(dz/100,1,d4s,d1s,d2s),2); %create dz grid the size of the NPP grid and convert units from cm to m
NPP_diat = cumsum(NPPv_diat.*dz_grid, 3); %mol C m-2 yr-1
NPP_diaz = cumsum(NPPv_diaz.*dz_grid, 3); %mol C m-2 yr-1
NPP_sp = cumsum(NPPv_sp.*dz_grid, 3); %mol C m-2 yr-1
NPP_sum = NPP_diat + NPP_diaz + NPP_sp;

% Air-sea fluxes
F_O2 = FG_O2.*secinyr/100/1000; %mmol O2/m^3 cm/sec to mol O2 m-2 yr-1
F_CO2 = FG_CO2.*secinyr/100/1000; %mmol CO2/m^3 cm/sec to mol CO2 m-2 yr-1

%% Interpolate data onto a finer depth grid
zmin = 5; zmax = 1000; zstep = 5; %define interpolated depth grid (in meters)
z_new = [zmin: zstep: zmax];

[POCflux_interp] = depth_interp_4D(z_top, z_new, POCflux);
[POCprod_interp] = depth_interp_4D(z, z_new, POCprod);
[POCremin_interp] = depth_interp_4D(z, z_new, POCremin);

[NPP_diat_interp] = depth_interp_4D(z, z_new, NPP_diat);
[NPP_diaz_interp] = depth_interp_4D(z, z_new, NPP_diaz);
[NPP_sp_interp] = depth_interp_4D(z, z_new, NPP_sp);
[NPP_sum_interp] = depth_interp_4D(z, z_new, NPP_sum);

%% Calculate depth criteria
%Euphotic depth, as defined by Lima et al. 2014: depth where POC production is equal to 1% of maximum POC production in the water column
[zeu] = zeuFromPOCProd(POCprod_interp,z_new);

%Particle compensation depth: depth of maximum POC flux
[zcomp, POCflux_zcomp] = zcompFromPOCFlux(z, POCflux);

%Calculate maximum annual MLD from seasonal MLD
mldmaxyr = NaN*ones(d1s,d2s,length(yrslist)); %initialize array with to hold mldmax for each year
mldmax = NaN*ones(d1s,d2s,d4s); %initialize array with time same size as input
for i = 1:length(yrslist)
    mldmaxyr(:,:,i) = max(mld(:,:,(i-1)*12+1:i*12),[],3);
    mldmax(:,:,(i-1)*12+1:i*12) = repmat(mldmaxyr(:,:,i),1,1,12);
end

%% Calculate values at depth criteria depths
%Create a cell array with all the depth criteria. 
    depthCriteria{1} = mld;
    depthCriteria{2} = zeu; 
    depthCriteria{3} = zcomp; 
    depthCriteria{4} = mldmax;
    depthCriteria{5} = 100*ones(d1s,d2s,d4s);
POCflux_depthCriteria = depthCriteriaVals(POCflux_interp, z_new, depthCriteria);
NPPsum_depthCriteria = depthCriteriaVals(NPP_sum_interp, z_new, depthCriteria);

%% Calculate climatological (over 10 years of output) and annual mean values

numCriteria = length(depthCriteria);
for i = 1:numCriteria
    %%%% For depth criteria
    depthCriteria_clim{i} = climatologicalMean(depthCriteria{i}); % Climatological monthly values (mean over the 10 years in dataset)
    depthCriteria_wAnnMean{i} = weightedAnnualMean(depthCriteria{i},... % Annual mean values for each year, weighted by monthly POC flux
        POCflux_depthCriteria{i}, year); 
    depthCriteria_wAnnMean_clim{i} = mean(depthCriteria_wAnnMean{i},3); % Climatological annual POC flux-weighted mean for each grid cell
    %%%% For POC flux rates
    POCflux_depthCriteria_clim{i} = climatologicalMean(POCflux_depthCriteria{i}); % Climatological monthly values (mean over the 10 years in dataset)
    for k = 1:length(yrslist) % Annual mean values for each year
        idyr = find(year == yrslist(k));
        POCflux_depthCriteria_AnnMean{i}(:,:,k) = mean(POCflux_depthCriteria{i}(:,:,idyr),3);
    end
    POCflux_depthCriteria_AnnMean_clim{i} = mean(POCflux_depthCriteria_AnnMean{i},3); % Climatological annual POC flux rates for each grid cell at each criterion
    %%%% For NPP rates
    NPPsum_depthCriteria_clim{i} = climatologicalMean(NPPsum_depthCriteria{i}); % Climatological monthly values (mean over the 10 years in dataset)
    for k = 1:length(yrslist) % Annual mean values for each year
        idyr = find(year == yrslist(k));
        NPPsum_depthCriteria_AnnMean{i}(:,:,k) = mean(NPPsum_depthCriteria{i}(:,:,idyr),3);
    end
    NPPsum_depthCriteria_AnnMean_clim{i} = mean(NPPsum_depthCriteria_AnnMean{i},3); % Climatological annual NPP rates for each grid cell at each criterion
end

%%%% For POC flux at all depths
POCflux_interp_clim = climatologicalMean(POCflux_interp);
for k = 1:length(yrslist)
    idyr = find(year == yrslist(k));
    POCflux_interp_AnnMean(:,:,:,k) = mean(POCflux_interp(:,:,:,idyr),4);
end
POCflux_interp_AnnMean_clim = mean(POCflux_interp_AnnMean,4);

%%%% For integrated NPP at all depths
NPPsum_interp_clim = climatologicalMean(NPP_sum_interp);
for k = 1:length(yrslist)
    idyr = find(year == yrslist(k));
    NPPsum_interp_AnnMean(:,:,:,k) = mean(NPP_sum_interp(:,:,:,idyr),4);
end
NPPsum_interp_AnnMean_clim = mean(NPPsum_interp_AnnMean,4);

%%%% For air-sea fluxes
F_O2_clim = climatologicalMean(F_O2);
F_CO2_clim = climatologicalMean(F_CO2);
for k = 1:length(yrslist)
    idyr = find(year == yrslist(k));
    F_O2_AnnMean(:,:,k) = mean(F_O2(:,:,idyr),3);
    F_CO2_AnnMean(:,:,k) = mean(F_CO2(:,:,idyr),3);
end
F_O2_AnnMean_clim = mean(F_O2_AnnMean,3);
F_CO2_AnnMean_clim = mean(F_CO2_AnnMean,3);

%% Calculate regional and global rates
    latint = 4;
latdiv = [-66:latint:66]; %define latitude divisions for zonal regions

clear NPP_sum_glob NPP_sum_reg POCflux_sum_glob POCflux_sum_reg area_reg
for i = 1:numCriteria
    [NPP_sum_glob(i), NPP_sum_reg(:,i), area_reg] = areaRateSum(NPPsum_depthCriteria_AnnMean_clim{i},latdiv,TLAT,TAREA,REGION_MASK);
    [POCflux_sum_glob(i), POCflux_sum_reg(:,i), area_reg] = areaRateSum(POCflux_depthCriteria_AnnMean_clim{i},latdiv,TLAT,TAREA,REGION_MASK);
end

[F_O2_sum_glob, F_O2_sum_reg, ~] = areaRateSum(F_O2_AnnMean_clim,latdiv,TLAT,TAREA,REGION_MASK);
[F_CO2_sum_glob, F_CO2_sum_reg, ~] = areaRateSum(F_CO2_AnnMean_clim,latdiv,TLAT,TAREA,REGION_MASK);

%% Calculate regional values for larger zones

% latdiv_broad = [-70,-30,-10,10,30,70]; %define latitude divisions for zonal regions
% 
% %clear NPP_sum_glob NPP_sum_reg POCflux_sum_glob POCflux_sum_reg area_reg
% for i = 1:numCriteria
%     [~, NPP_sum_reg_broad(:,i), area_reg] = areaRateSum(NPPsum_depthCriteria_AnnMean_clim{i},latdiv_broad,TLAT,TAREA,REGION_MASK);
%     [~, POCflux_sum_reg_broad(:,i), area_reg] = areaRateSum(POCflux_depthCriteria_AnnMean_clim{i},latdiv_broad,TLAT,TAREA,REGION_MASK);
% end
% 
% POCflux_highlat = sum(POCflux_sum_reg_broad([1,5],:).*repmat((area_reg([1,5]))/(area_reg(1)+area_reg(5)),5,1)'); %weighted mean by area
% POCflux_subtrop = sum(POCflux_sum_reg_broad([2,4],:).*repmat((area_reg([2,4]))/(area_reg(2)+area_reg(4)),5,1)'); %weighted mean by area
% NPP_highlat = sum(NPP_sum_reg_broad([1,5],:).*repmat((area_reg([1,5]))/(area_reg(1)+area_reg(5)),5,1)'); %weighted mean by area
% NPP_subtrop = sum(NPP_sum_reg_broad([2,4],:).*repmat((area_reg([2,4]))/(area_reg(2)+area_reg(4)),5,1)'); %weighted mean by area
% 
% eratio_highvslowlat = (POCflux_highlat./NPP_highlat)./(POCflux_subtrop./NPP_subtrop);

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
lon_gridsize = 4; lat_gridsize = 4;

[glon, glat, POCflux_interp_AnnMean_clim_grid] = regrid_even(TLONG, TLAT, POCflux_interp_AnnMean_clim,...
    lonmin, lonmax, latmin, latmax, lon_gridsize, lat_gridsize);
[glon, glat, NPPsum_interp_AnnMean_clim_grid] = regrid_even(TLONG, TLAT, NPPsum_interp_AnnMean_clim,...
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
stnNames = {'Drake Passage (north)','Irminger Sea','HOT','OSP','Kuroshio','Equatorial Pacific','W. Subtropical N. Pacific','Subtropical S. Atlantic','Southern Ocean'};
zmaxstn = [350 850 150 150 350 150 150 150 450];

%% Run plotting script to make figures
CSMBGC_plotting

%% Calculate rate of attenuation below compensation depth
%Calculating just for climatological annual mean
ZComp = depthCriteria_wAnnMean_clim{3};
depthCriteria_belowZComp{1} = ZComp - 50;
depthCriteria_belowZComp{2} = ZComp - 100;
POCflux_belowZComp = depthCriteriaVals(POCflux_interp_AnnMean_clim, z_new, depthCriteria_belowZComp);

figure(100); clf
for i = 1:2
subplot(2,1,i)
    Remin_fract = 1 - (POCflux_belowZComp{i}./POCflux_depthCriteria_AnnMean_clim{3}); %Fraction remineralized in 1st X m below ZComp
    histogram(Remin_fract(:))
    title([num2str(nanmean(Remin_fract(:))*100) '% +/- ' num2str(nanstd(Remin_fract(:))*100) ' remineralized'])
end

%% Extract station data to use as forcing/validation data for box model

% Use nearest point in CSM grid for each of the time series stations
% (inds = coordinates in grid; gridpt = lat/lon values of grid pts)
for i = 1:5
    SST_stn = squeeze(TEMP(inds(i,1),inds(i,2),1,:)); %deg C
    SSS_stn = squeeze(SALT(inds(i,1),inds(i,2),1,:)); %psu
    DICsurf_stn = squeeze(DIC(inds(i,1),inds(i,2),1,:)); %mmol C/m^3
    ALKsurf_stn = squeeze(ALK(inds(i,1),inds(i,2),1,:)); %meq C/m^3
    mld_stn = squeeze(mld(inds(i,1),inds(i,2),:)); %m
    PV_CO2_stn = squeeze(PV_CO2(inds(i,1),inds(i,2),:))/100*(60*60*24); %cm/sec --> m/day
    pCO2surf_stn = squeeze(pCO2SURF(inds(i,1),inds(i,2),:)); %ppmv
    FG_CO2_stn = squeeze(FG_CO2(inds(i,1),inds(i,2),:))/100*(60*60*24); %mmol DIC/m^3 cm/sec --> mmol DIC/m^2/day 
    POCflux_fromML_stn = squeeze(POCflux_depthCriteria{1}(inds(i,1),inds(i,2),:)); %mol C m-2 yr-1
    POCflux_MLDmax_stn = squeeze(POCflux_depthCriteria{4}(inds(i,1),inds(i,2),:)); %mol C m-2 yr-1
    POCremin_stn = squeeze(POCremin_interp(inds(i,1),inds(i,2),:,:)); %keep full depth (depths are z_new) - mol C m-3 yr-1
    POCprod_stn = squeeze(POCprod_interp(inds(i,1),inds(i,2),:,:)); %keep full depth (depths are z_new) - mol C m-3 yr-1


%Calculate flux of remineralized POC into mixed layer from material
%entrained from thermocline
remin_therm = zeros(length(z_new),length(time));
remin_entrain = zeros(length(time),1);
for j = 2:length(time)
    indtherm = find(z_new < max(mld_stn) & z_new >= mld_stn(j)); %find depths within seasonal thermocline
        remin_therm(indtherm,j) = remin_therm(indtherm,j-1) + POCremin_stn(indtherm,j) - POCprod_stn(indtherm,j); %accumulate remineralized POC in thermocline (remineralization - production)
    if mld_stn(j) > mld_stn(j-1)
        indent = find(z_new > mld_stn(j-1) & z_new <= mld_stn(j)); %find depths just entrained into mixed layer
        remin_entrain(j) = sum(remin_therm(indent,j-1)*(z_new(2) - z_new(1))); %integrated entrainment through column just entrained - mol C m-2 yr-1
        remin_therm(indent,j) = 0;
    end
end

figure(200 + i); clf
subplot(4,2,1); plot(SST_stn); title(['SST (deg C) at ' stnNames{i}])
subplot(4,2,2); plot(SSS_stn); title(['SSS (psu) at ' stnNames{i}])
subplot(4,2,3); plot(DICsurf_stn); title(['DIC (mmol m^{-3}) at ' stnNames{i}])
subplot(4,2,4); plot(pCO2surf_stn); title(['pCO_2 (ppmv) at ' stnNames{i}])
subplot(4,2,5); plot(mld_stn); set(gca,'ydir','reverse'); title(['MLD (m) at ' stnNames{i}])
subplot(4,2,6); plot(PV_CO2_stn); title(['Piston velocity for CO_2 (m d^{-1}) at ' stnNames{i}])
subplot(4,2,7); plot(-FG_CO2_stn);  title(['Air-sea CO_2 flux (mmol m^{-2} yr^{-1}) at ' stnNames{i}])
subplot(4,2,8); plot(POCflux_fromML_stn); hold on; plot(-remin_entrain,'c'); plot(POCflux_fromML_stn - remin_entrain,'k'); hold on; title(['POC flux (from ML, back to ML, and net; mmol m^{-2} yr^{-1}) at ' stnNames{i}])

    stn_location = gridpt(i,:);
save(['CCSMsurfdata_' stnNames{i}], 'stn_location', 'SST_stn', 'SSS_stn', 'DICsurf_stn', 'ALKsurf_stn', 'pCO2surf_stn', 'mld_stn',...
    'PV_CO2_stn', 'FG_CO2_stn', 'POCflux_fromML_stn', 'remin_entrain')

end

