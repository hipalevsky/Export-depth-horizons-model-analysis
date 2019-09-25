%% Make surface plots and sections of CSMBGC output
%First run CSMBGC_openfiles

%Make all of these rough plots on the CSM grid --> later will want to
%regrid everything

%%%% In comparing POC flux at different depth criteria, key to consider will be:
% 1) how much POC_PROD was missed using a given depth?
% 2) how much POC_REMIN happened below a given depth but above others?

%% Convert times
    year = floor((time-1)/365); %the "-1" keeps last day of Dec in current yr
        yrslist = unique(year);
    day = (time-1) - year*365 + 1;
    
%% Calculate POC flux from winter MLD, euphotic depth (defined in Lima 2014
% as depth where POC production is equal to 1% of maximum POC production in the water column), 100 m, and seasonal MLD and make maps
    z = z_t/100; %depth from surface to midpoint of layer in meters
    mld = HMXL/100; %mixed layer depth in meters

%Calculate maximum annual mixed layer depth and POC flux at maximum annual MLD
[m,n] = size(TLAT); %size of spatial grid    
mldmax = NaN*ones(m,n,length(yrslist));
POCflux_mldmax = mldmax;
z_id_mldmax = mldmax; %find nearest depth interval to mld *****(should consider other options for interpolating within depth grids)
for i = 1:length(yrslist)
    idyr = find(year == yrslist(i));
    for j = 1:m
        for k = 1:n
            mldmax(j,k,i) = max(mld(j,k,idyr)); %maximum annual MLD
            if isnan(mldmax(j,k,i)) == 0
                [~,z_id_mldmax(j,k,i)] = min(abs(z - mldmax(j,k,i))); %nearest depth (defined at middle of grid) to maximum annual MLD
                POCflux_mldmax(j,k,i) = mean(POC_FLUX_IN(j,k,z_id_mldmax(j,k,i),idyr))*60*60*24/100; %annual mean POC flux at maximum annual MLD in mmol C m-2 d-1
            end
        end
    end
end
%%
%Calculate euphotic depth (as POC production is equal to 1% of maximum POC
%production in the water column) and POC flux at mean annual euphotic depth
z_eu = NaN*ones(m,n,length(time));
POCflux_z_eu = z_eu;
POCflux_mld = z_eu;
z_comp = z_eu;
POCflux_comp = z_eu;
for i = 1:length(time)
    for j = 1:m
        for k = 1:n
            %POC flux at variable euphotic depth
            POC_max = max(POC_PROD(j,k,:,i)); %maximum POC production in water column for given location/month
            if isnan(POC_max) == 0
                [~,z_eu_id(j,k,i)] = min(abs(POC_PROD(j,k,:,i) - POC_max/100)); %nearest depth (defined at middle of grid) to 1% of maximum POC production
                z_eu(j,k,i) = z(z_eu_id(j,k,i));
                POCflux_z_eu(j,k,i) = POC_FLUX_IN(j,k,z_eu_id(j,k,i),i)*60*60*24/100; %POC flux at seasonal z_eu in mmol C m-2 d-1
            end
            %POC flux at variable seasonal MLD
            if isnan(mld(j,k,i)) == 0
                [~,mld_id(j,k,i)] = min(abs(z - mld(j,k,i))); %nearest depth (defined at middle of grid) to seasonal MLD
                POCflux_mld(j,k,i) = POC_FLUX_IN(j,k,mld_id(j,k,i),i)*60*60*24/100; %POC flux at seasonal mld in mmol C m-2 d-1
            end
            %Calculate POC flux at particle compensation depth
            [POC_flux_max,z_comp_id] = max(POC_FLUX_IN(j,k,:,i)); %maximum POC flux in water column for given location/month
            if isnan(POC_flux_max) == 0
                z_comp(j,k,i) = z(z_comp_id); %particle compensation depth
                POCflux_comp(j,k,i) = POC_flux_max*60*60*24/100; %POC flux at particle compensation depth in mmol C m-2 d-1
            end
        end
    end
end

%% Alternative calculation of euphotic zone depth and export flux (using annual mean over all 10 years)
POC_PROD_mean = mean(POC_PROD,4); %mean POC production over all months of all 10 years
z_eu_mean = NaN*(ones(m,n)); z_eu_mean_id = NaN*(ones(m,n)); POCflux_z_eu_mean = NaN*(ones(m,n));
for j = 1:m
    for k = 1:n
        POC_max = max(POC_PROD_mean(j,k,:));
        if isnan(POC_max) == 0
            [~,z_eu_mean_id(j,k)] = min(abs(POC_PROD_mean(j,k,:) - POC_max/100)); %nearest depth (defined at middle of grid) to 1% of maximum POC production
            z_eu_mean(j,k) = z(z_eu_mean_id(j,k));
            if isnan(z_eu_mean(j,k)) == 0
                POCflux_z_eu_mean(j,k) = mean(POC_FLUX_IN(j,k,z_eu_mean_id(j,k),:))*60*60*24/100; %annual mean POC flux at maximum annual MLD in mmol C m-2 d-1
            end
        end
    end
end

%% Calculate annual POC fluxes at the base of the seasonally varying mixed layer depth/euphotic depth/100 m

for i = 1:length(yrslist)
    idyr = find(year == yrslist(i));
    POCflux_100m_ann(:,:,i) = mean(squeeze(POC_FLUX_IN(:,:,9,idyr))*60*60*24/100,3); %POC flux at 100 m in mmol C m-2 d-1
    POCflux_z_eu_ann(:,:,i) = mean(squeeze(POCflux_z_eu(:,:,idyr)),3);
    POCflux_seasmld_ann(:,:,i) = mean(squeeze(POCflux_mld(:,:,idyr)),3);
    POCflux_comp_ann(:,:,i) = mean(squeeze(POCflux_comp(:,:,idyr)),3);
end

%% Calculate mean annual euphotic, compensation, and mixed layer depths, weighted by POC flux
z_eu_wmean = NaN*(ones(m,n,length(yrslist))); z_comp_wmean = NaN*(ones(m,n,length(yrslist))); seasmld_wmean = NaN*(ones(m,n,length(yrslist))); 

for j = 1:m
    for k = 1:n
        for i = 1:length(yrslist)
            idyr = find(year == yrslist(i));
                fract_eu = squeeze(POCflux_z_eu(j,k,idyr))./(12*POCflux_z_eu_ann(j,k,i));
            z_eu_wmean(j,k,i) = sum(squeeze(z_eu(j,k,idyr)).*fract_eu);
                fract_comp = squeeze(POCflux_comp(j,k,idyr))./(12*POCflux_comp_ann(j,k,i));
            z_comp_wmean(j,k,i) = sum(squeeze(z_comp(j,k,idyr)).*fract_comp);
                fract_mld = squeeze(POCflux_mld(j,k,idyr))./(12*POCflux_seasmld_ann(j,k,i));
            seasmld_wmean(j,k,i) = sum(squeeze(mld(j,k,idyr)).*fract_mld);
            mldmin(j,k,i) = min(squeeze(mld(j,k,idyr)));
        end
    end
end

%% Calculate POC flux using each depth criterion (global/regional)
    moltoPg = 12*10^-15;

AA = (mean(POCflux_100m_ann,3)*365/1000).*TAREA/(100^2)*moltoPg; %POC flux in each cell (mol C m-2 yr-1) * area of each cell (converted from cm^2 to m^2) --> mol C yr-1 in each cell
BB = (mean(POCflux_seasmld_ann,3)*365/1000).*TAREA/(100^2)*moltoPg; %POC flux in each cell (mol C m-2 yr-1) * area of each cell (converted from cm^2 to m^2) --> mol C yr-1 in each cell
CC = ((POCflux_z_eu_mean)*365/1000).*TAREA/(100^2)*moltoPg; %POC flux in each cell (mol C m-2 yr-1) * area of each cell (converted from cm^2 to m^2) --> mol C yr-1 in each cell
DD = (mean(POCflux_mldmax,3)*365/1000).*TAREA/(100^2)*moltoPg; %POC flux in each cell (mol C m-2 yr-1) * area of each cell (converted from cm^2 to m^2) --> mol C yr-1 in each cell
EE = (mean(POCflux_comp_ann,3)*365/1000).*TAREA/(100^2)*moltoPg; %POC flux in each cell (mol C m-2 yr-1) * area of each cell (converted from cm^2 to m^2) --> mol C yr-1 in each cell
FF = (mean(POCflux_z_eu_ann,3)*365/1000).*TAREA/(100^2)*moltoPg; %POC flux in each cell (mol C m-2 yr-1) * area of each cell (converted from cm^2 to m^2) --> mol C yr-1 in each cell

    reglist = unique(REGION_MASK); %find all region mask numbers
%Initialize arrays to hold regional and global sums and calculate
%global values
    POCflux_100m_sums = NaN*ones(length(reglist) + 1,1);
        POCflux_100m_sums(1) = nansum(nansum(AA));
    POCflux_seasmld_sums = NaN*ones(length(reglist) + 1,1);
        POCflux_seasmld_sums(1) = nansum(nansum(BB));
    POCflux_zeumean_sums = NaN*ones(length(reglist) + 1,1);
        POCflux_zeumean_sums(1) = nansum(nansum(CC));
    POCflux_mldmax_sums = NaN*ones(length(reglist) + 1,1);
        POCflux_mldmax_sums(1) = nansum(nansum(DD));
    POCflux_comp_sums = NaN*ones(length(reglist) + 1,1);
        POCflux_comp_sums(1) = nansum(nansum(EE));
    POCflux_zeuseas_sums = NaN*ones(length(reglist) + 1,1);
        POCflux_zeuseas_sums(1) = nansum(nansum(FF));
%Calculate values for each region
for i = 1:length(reglist);
    for j = 1:m
        idreg = find(REGION_MASK(j,:) == reglist(i));
        aa(j) = nansum(AA(j,idreg));
        bb(j) = nansum(BB(j,idreg));
        cc(j) = nansum(CC(j,idreg));
        dd(j) = nansum(DD(j,idreg));
        ee(j) = nansum(EE(j,idreg));
        ff(j) = nansum(FF(j,idreg));
    end
    POCflux_100m_sums(i+1) = nansum(aa);
    POCflux_seasmld_sums(i+1) = nansum(bb);
    POCflux_zeumean_sums(i+1) = nansum(cc);
    POCflux_mldmax_sums(i+1) = nansum(dd);
    POCflux_comp_sums(i+1) = nansum(ee);
    POCflux_zeuseas_sums(i+1) = nansum(ff);
end
    

%% Plot of POC flux using each depth criterion
figure(4); clf
cmin = 0; cmax = 6;
C = flipud(colormap(lbmap(61,'RedBlue')));  C(1,:) = [0 0 0];
    subplot(3,2,1)
imagesc(rot90(mean(POCflux_100m_ann,3)*365/1000)); hc = colorbar; caxis([cmin cmax]); colormap(C);
title(['Annual POC flux at 100 m, ' num2str(POCflux_100m_sums(1),2) ' Pg C yr^{-1}']); ylabel(hc,'mol C m^{-2} yr^{-1}');
    subplot(3,2,2)
imagesc(rot90(mean(POCflux_seasmld_ann,3)*365/1000)); hc = colorbar; caxis([cmin cmax]); colormap(C);
title(['Annual POC flux at seasonally-varying MLD, ' num2str(POCflux_seasmld_sums(1),2) ' Pg C yr^{-1}']); ylabel(hc,'mol C m^{-2} yr^{-1}');
    subplot(3,2,3)
imagesc(rot90(POCflux_z_eu_mean*365/1000)); hc = colorbar; caxis([cmin cmax]); colormap(C);
title(['Annual POC flux at mean annual euphotic depth, ' num2str(POCflux_zeumean_sums(1),2) ' Pg C yr^{-1}']); ylabel(hc,'mol C m^{-2} yr^{-1}');
    subplot(3,2,4)
imagesc(rot90(mean(POCflux_mldmax,3)*365/1000)); hc = colorbar; caxis([cmin cmax]); colormap(C);
title(['Annual POC flux at Winter MLD, ' num2str(POCflux_mldmax_sums(1),2) ' Pg C yr^{-1}']); ylabel(hc,'mol C m^{-2} yr^{-1}');
    subplot(3,2,5)
imagesc(rot90(mean(POCflux_z_eu_ann,3)*365/1000)); hc = colorbar; caxis([cmin cmax]); colormap(C);
title(['Annual POC flux at seasonally-varying euphotic depth, ' num2str(POCflux_zeuseas_sums(1),2) ' Pg C yr^{-1}']); ylabel(hc,'mol C m^{-2} yr^{-1}');
    subplot(3,2,6)
imagesc(rot90(mean(POCflux_comp_ann,3)*365/1000)); hc = colorbar; caxis([cmin cmax]); colormap(C);
title(['Annual POC flux at seasonally-varying particle compensation depth, ' num2str(POCflux_comp_sums(1),2) ' Pg C yr^{-1}']); ylabel(hc,'mol C m^{-2} yr^{-1}');

%% Map of depth criteria
figure(40); clf
C = colormap(parula);  C(1,:) = [0 0 0];
cmin = 0; cmax = 300;
    subplot(2,1,1)
imagesc(rot90(mean(mldmax,3))); hc = colorbar; caxis([cmin cmax]); colormap(C); title('Winter MLD'); ylabel(hc,'Depth (m)');
    subplot(2,1,2)
imagesc(rot90(z_eu_mean)); hc = colorbar; caxis([cmin cmax]); colormap(C); title('Mean annual euphotic depth'); ylabel(hc,'Depth (m)');

figure(41); clf
C2 = colormap(lbmap(64,'RedBlue'));  C2(1,:) = [0 0 0];
cmin2 = -100; cmax2 = 100;
    imagesc(rot90(mean(mldmax,3)-z_eu_mean)); hc = colorbar; caxis([cmin2 cmax2]); colormap(C2); title('Winter MLD - Mean annual euphotic depth'); ylabel(hc,'Depth difference (m)');


%% POC flux differences among depth criteria
figure(5); clf
cmin = -1; cmax = 1;
C = colormap(lbmap(64,'RedBlue'));  C(1,:) = [0 0 0];
    subplot(3,2,1)
imagesc(rot90((mean(POCflux_seasmld_ann,3)-mean(POCflux_100m_ann,3))*-365/1000)); hc = colorbar; caxis([cmin cmax]); colormap(C); title('Annual POC flux: 100 m - Seasonally-varying MLD'); ylabel(hc,'mol C m^{-2} yr^{-1}');
    subplot(3,2,3)
imagesc(rot90((POCflux_z_eu_mean-mean(POCflux_100m_ann,3))*-365/1000)); hc = colorbar; caxis([cmin cmax]); colormap(C); title('Annual POC flux: 100 m - Mean annual euphotic depth'); ylabel(hc,'mol C m^{-2} yr^{-1}');
    subplot(3,2,5)
imagesc(rot90((mean(POCflux_mldmax,3)-mean(POCflux_100m_ann,3))*-365/1000)); hc = colorbar; caxis([cmin cmax]); colormap(C); title('Annual POC flux: 100 m - Winter MLD'); ylabel(hc,'mol C m^{-2} yr^{-1}');
    subplot(3,2,2)
imagesc(rot90((mean(POCflux_comp_ann,3)-mean(POCflux_100m_ann,3))*-365/1000)); hc = colorbar; caxis([cmin cmax]); colormap(C); title('Annual POC flux: 100 m - Seasonally-varying compensation depth'); ylabel(hc,'mol C m^{-2} yr^{-1}');
    subplot(3,2,4)
imagesc(rot90((mean(POCflux_z_eu_ann,3)-mean(POCflux_100m_ann,3))*-365/1000)); hc = colorbar; caxis([cmin cmax]); colormap(C); title('Annual POC flux: Seasonally-varying euphotic depth - Mean annual euphotic depth'); ylabel(hc,'mol C m^{-2} yr^{-1}');
    subplot(3,2,6)
imagesc(rot90((mean(POCflux_mldmax,3)-mean(POCflux_z_eu_ann,3))*-365/1000)); hc = colorbar; caxis([cmin cmax]); colormap(C); title('Annual POC flux: 100 m - Seasonally-varying euphotic depth'); ylabel(hc,'mol C m^{-2} yr^{-1}');

%% Depth dependence of differences among depth criteria
figure(50); clf
    A = -(mean(POCflux_mldmax,3)-mean(POCflux_100m_ann,3)); A = A(:);
    B = mean(mldmax,3); B = B(:);
    ind = find(B > 50 & isnan(A + B) == 0);
        idsp = [1:1:(length(ind))]; %subsample one in x points
    [rho_100,df_100,rho_sig95_100] = correlate(B(ind(idsp)),A(ind(idsp)));
dscatter(B(ind(idsp)),A(ind(idsp))); hold on;
plot([1:1000],0*ones(1000,1),'k--'); hold on;
axis([45 500 -4 8])
xlabel('Maximum annual MLD (m)')
ylabel('Annual POC flux offset (mol C m^{-2} yr^{-1})')
title(['POC flux, 100 m - Winter MLD, r^2 = ' num2str(rho_100^2,2)])

%% POC flux section
    figure(7); clf
POC_flux_mean = mean(POC_FLUX_IN,4)*60*60*24/100*365/1000; %mean POC flux over all months of all 10 years in mol C m-2 yr-1
%Set parameters for plotting
lonmin = 2; lonmax = 115;
zmin = 0; zmax = 36;
cmin = 0; cmax = 5;
L = 2;
Atl = 5; Pac = 68; %slices to use for sections
C = (colormap(lbmap(61,'Blue'))); C(1,:) = [1, 1, 1];

%Interpolate data for sections onto an even depth grid
    zgrid = [0:10:500];
    Atl_sect = NaN*ones(116,length(zgrid)); Pac_sect = NaN*ones(116,length(zgrid));
    for i = 1:116
        Atl_sect(i,:) = interp1(z,squeeze(POC_flux_mean(Atl,i,:)),zgrid,'linear');
        Pac_sect(i,:) = interp1(z,squeeze(POC_flux_mean(Pac,i,:)),zgrid,'linear');
    end
%Calculate MLD_max and mean z_eu for sections
mldmax_clim = mean(mldmax,3);
z_eu_clim = mean(z_eu,3);
%mldmin = min(mld,[],3);
    z_eu_min = min(z_eu,[],3);
    z_eu_max = max(z_eu,[],3);
    z_comp_min = min(z_comp,[],3);
    z_comp_max = max(z_comp,[],3);
    
    subplot(2,1,1)
imagesc(flipud(rot90(Atl_sect))); hc = colorbar; caxis([cmin cmax]); colormap(C); hold on;
plot([1:116],11*ones(116,1),'w','linewidth',L); hold on;
plot([1:116],mldmax_clim(Atl,:)/10+1,'k','linewidth',L); hold on;
%plot([1:116],z_eu_mean(Atl,:)/10+1,'color',nicecolor('ryw'),'linewidth',L); hold on;
plot([1:116],mean(z_eu_wmean(Atl,:,:),3)/10+1,'color',nicecolor('ryw'),'linewidth',L); hold on;
    %plot([1:116],z_eu_min(Atl,:)/10+1,'--','color',nicecolor('ryw'),'linewidth',L); hold on;
    %plot([1:116],z_eu_max(Atl,:)/10+1,'--','color',nicecolor('ryw'),'linewidth',L); hold on;
plot([1:116],mean(mldmin(Atl,:,:),3)/10+1,'color',nicecolor('m'),'linewidth',L); hold on;
    %plot([1:116],mean(seasmld_wmean(Atl,:,:),3)/10+1,'--','color',nicecolor('m'),'linewidth',L); hold on;
plot([1:116],mean(z_comp_wmean(Atl,:,:),3)/10+1,'color',nicecolor('ccgk'),'linewidth',L); hold on;
    %plot([1:116],z_comp_min(Atl,:)/10+1,'--','color',nicecolor('c'),'linewidth',L); hold on;
    %plot([1:116],z_comp_max(Atl,:)/10+1,'--','color',nicecolor('c'),'linewidth',L); hold on;
title('Atlantic section (5), South <--> North'); ylabel(hc,'mol C m^{-2} yr^{-1}'); set(gca,'Ytick',[11:10:51],'YTickLabel',{'100','200','300','400','500'}); ylabel('Depth (m)');
axis([lonmin lonmax zmin zmax])
    subplot(2,1,2)
imagesc(flipud(rot90(Pac_sect))); hc = colorbar; caxis([cmin cmax]); colormap(C); hold on;
plot([1:116],11*ones(116,1),'w','linewidth',L); hold on;
plot([1:116],mldmax_clim(Pac,:)/10+1,'k','linewidth',L); hold on;
%plot([1:116],z_eu_mean(Pac,:)/10+1,'--','color',nicecolor('ryw'),'linewidth',L); hold on;
plot([1:116],mean(z_eu_wmean(Pac,:,:),3)/10+1,'color',nicecolor('ryw'),'linewidth',L); hold on;
    %plot([1:116],z_eu_min(Pac,:)/10+1,'--','color',nicecolor('ryw'),'linewidth',L); hold on;
    %plot([1:116],z_eu_max(Pac,:)/10+1,'--','color',nicecolor('ryw'),'linewidth',L); hold on;
plot([1:116],mean(mldmin(Pac,:,:),3)/10+1,'color',nicecolor('m'),'linewidth',L); hold on;
    %plot([1:116],mean(seasmld_wmean(Pac,:,:),3)/10+1,'--','color',nicecolor('m'),'linewidth',L); hold on;
plot([1:116],mean(z_comp_wmean(Pac,:,:),3)/10+1,'color',nicecolor('ccgk'),'linewidth',L); hold on;
    %plot([1:116],z_comp_min(Pac,:)/10+1,'--','color',nicecolor('c'),'linewidth',L); hold on;
    %plot([1:116],z_comp_max(Pac,:)/10+1,'--','color',nicecolor('c'),'linewidth',L); hold on;
title('Pacific section (68), South <--> North'); ylabel(hc,'mol C m^{-2} yr^{-1}'); set(gca,'Ytick',[11:10:51],'YTickLabel',{'100','200','300','400','500'}); ylabel('Depth (m)');
axis([lonmin lonmax zmin zmax])

%% Look systematically at relationship among all depth criteria
%Vizualize the comparison: mimimum annual MLD < compensation depth <
%euphotic depth in nearly all cases
figure(10); clf
C = colormap(lbmap(61,'RedBlue')); C(1,:) = [0 0 0];
cmin = -30; cmax = 30;
    subplot(311)
imagesc(rot90(mean(z_eu_wmean,3)-mean(z_comp_wmean,3))); colorbar; colormap(C); caxis([cmin cmax]);
title('Euphotic depth - Compensation depth (m)');
    subplot(312)
imagesc(rot90(mean(z_comp_wmean,3)-mean(mldmin,3))); colorbar; colormap(C); caxis([cmin cmax]);
title('Compensation depth - Minimum annual MLD (m)');
    subplot(313)
imagesc(rot90(mean(mldmin,3))); colorbar; colormap(C); caxis([0 100]);
title('Minimum annual MLD (m)');

%% Calculate relationships with maximum annual MLD to create regions
depthregID = zeros(m,n,length(yrslist));
for j = 1:m
    for k = 1:n
        for i = 1:length(yrslist)
            %Assign regions based on relationship between maximum annual
            %MLD and other depth criteria
            if mldmax(j,k,i) > z_eu_wmean(j,k,i)
                depthregID(j,k,i) = 1; %maximum mld > euphotic depth
            elseif mldmax(j,k,i) <= z_eu_wmean(j,k,i) & mldmax(j,k,i) > z_comp_wmean(j,k,i);
                depthregID(j,k,i) = 2; %euphotic depth > maximum mld > compensation depth
            elseif mldmax(j,k,i) < z_comp_wmean(j,k,i)
                depthregID(j,k,i) = 3;
            end
            %Note locations where compensation depth < minimum annual MLD
            if z_comp_wmean(j,k,i) < mldmin(j,k,i)
                depthregID(j,k,i) = depthregID(j,k,i) + 0.5;
            end
        end
    end
end

% Plot initial regions
figure(11); clf
    C = [nicecolor('k'); nicecolor('ry');  nicecolor('ggk'); nicecolor('bbc'); nicecolor('rb')];
    imagesc(rot90([mean(depthregID,3); mean(depthregID,3)])); colorbar; colormap(C);

% Assign all points where compensation depth < minimum annual MLD or where
% there is interannual variation in region assignment to an "unassigned"
% region
depthregID_clim = mean(depthregID,3);
for i = 1:m
    for j = 1:n
        if depthregID_clim(i,j) - floor(depthregID_clim(i,j)) > 0
            depthregID_clim(i,j) = 4;
        end
    end
end
depthregID_clim_double = [depthregID_clim; depthregID_clim];

% Plot climatologial regions
figure(12); clf
    C = [nicecolor('k'); nicecolor('ry');  nicecolor('bbc'); nicecolor('rb'); nicecolor('wwk')];
    imagesc(rot90(depthregID_clim_double(45:155,:))); hc = colorbar('southoutside'); colormap(C);
    set(gca,'Xtick',[],'Ytick',[]);
    set(hc,'Xtick',[1 2 3],'Xticklabel',{'MLDmax > zeu','zeu > MLDmax > zcomp','MLDmax < zcomp'},'Fontsize',10);
%% Calculate POC flux using each depth criterion (global/regional)
    moltoPg = 12*10^-15;

AA = (mean(POCflux_100m_ann,3)*365/1000).*TAREA/(100^2)*moltoPg; %POC flux in each cell (mol C m-2 yr-1) * area of each cell (converted from cm^2 to m^2) --> mol C yr-1 in each cell
BB = (mean(POCflux_seasmld_ann,3)*365/1000).*TAREA/(100^2)*moltoPg; %POC flux in each cell (mol C m-2 yr-1) * area of each cell (converted from cm^2 to m^2) --> mol C yr-1 in each cell
CC = ((POCflux_z_eu_mean)*365/1000).*TAREA/(100^2)*moltoPg; %POC flux in each cell (mol C m-2 yr-1) * area of each cell (converted from cm^2 to m^2) --> mol C yr-1 in each cell
DD = (mean(POCflux_mldmax,3)*365/1000).*TAREA/(100^2)*moltoPg; %POC flux in each cell (mol C m-2 yr-1) * area of each cell (converted from cm^2 to m^2) --> mol C yr-1 in each cell
EE = (mean(POCflux_comp_ann,3)*365/1000).*TAREA/(100^2)*moltoPg; %POC flux in each cell (mol C m-2 yr-1) * area of each cell (converted from cm^2 to m^2) --> mol C yr-1 in each cell
FF = (mean(POCflux_z_eu_ann,3)*365/1000).*TAREA/(100^2)*moltoPg; %POC flux in each cell (mol C m-2 yr-1) * area of each cell (converted from cm^2 to m^2) --> mol C yr-1 in each cell

    reglist = [1,2,3]; %three depth-based regions
%Initialize arrays to hold regional and global sums and calculate global values
    POCfluxsum_depthreg = NaN*ones(length(reglist),6);
%Calculate values for each region
for i = 1:length(reglist);
    for j = 1:m
        idreg = find(depthregID_clim(j,:) == reglist(i));
        aa(j) = nansum(AA(j,idreg));
        bb(j) = nansum(BB(j,idreg));
        cc(j) = nansum(CC(j,idreg));
        dd(j) = nansum(DD(j,idreg));
        ee(j) = nansum(EE(j,idreg));
        ff(j) = nansum(FF(j,idreg));
        area(j) = nansum(TAREA(j,idreg));
    end
    POCfluxsum_depthreg(i,1) = nansum(aa);
    POCfluxsum_depthreg(i,2) = nansum(bb);
    POCfluxsum_depthreg(i,3) = nansum(cc);
    POCfluxsum_depthreg(i,4) = nansum(dd);
    POCfluxsum_depthreg(i,5) = nansum(ee);
    POCfluxsum_depthreg(i,6) = nansum(ff);
    TAREA_depthreg(i) = nansum(area);
end

% Calculate percent of POC flux from each region with each depth criterion
POCfluxtotal = [POCflux_100m_sums(1) POCflux_seasmld_sums(1) POCflux_zeumean_sums(1)...
    POCflux_mldmax_sums(1) POCflux_comp_sums(1) POCflux_zeuseas_sums(1)];
POCfluxpercent_depthreg = POCfluxsum_depthreg./repmat(POCfluxtotal,3,1)*100;
TAREApercent_depthreg = TAREA_depthreg./sum(TAREA_depthreg);

%% Pull out data and plot points of time series stations

%HOT, BATS, Irminger Sea
stnlat = [22.75 31.67 60];
stnlon = [-158+360 360-(64+10/60) 360-39.5];

tol = 3;
for i = 1
    latdif = abs(TLAT - stnlat(i));
    londif = abs(TLONG - stnlon(i));
    indlat = find(latdif < tol);
    indlon = find(londif < tol);
end


%NEXT STEP WILL THEN BE TO LOOK AT CHARACTERISTICS OF THE THREE REGION
%TYPES, AND TO PULL OUT INDIVIDUAL "STATIONS" WITHIN EACH REGION

