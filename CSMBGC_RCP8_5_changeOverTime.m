% Script to analyze output from different time periods within 21st century,
% extracted from RCP 8.5 scenario output in CSMBGC_exportanalysis_RCP8_5

%% Add external hard drive to path
addpath('D:/b40.rcp8_5.1deg.bdrd.001')  
addpath('D:/b40.rcp8_5.1deg.bdrd.001.POC')  

%% Add previously processed output to path
addpath('C:/Users/palevsky/Dropbox/MATLAB/CSMBGC data processing/RCP 8_5 processed output')

%% Load constantly relevant variables by loading data from a single year
[z, z_top, TLAT, TLONG, TAREA, REGION_MASK, ~, ~, ~, ~, ~, ~, ~] = CSMBGC_openfiles_RCP8_5(2100);

%% Load output from both beginning and end of century, then wrap each into a structure
load RCP8_5_2005to2024_ClimOnly.mat

outputToday = v2struct(yrslist, depthCriteria_clim, depthCriteria_wAnnMean, depthCriteria_wAnnMean_clim,...
    POCflux_depthCriteria_clim, POCflux_depthCriteria_AnnMean, POCflux_depthCriteria_AnnMean_clim,...
    NPPsum_depthCriteria_clim, NPPdiat_depthCriteria_clim, NPPdiaz_depthCriteria_clim, NPPsp_depthCriteria_clim,...
    NPPsum_depthCriteria_AnnMean, NPPdiat_depthCriteria_AnnMean, NPPdiaz_depthCriteria_AnnMean, NPPsp_depthCriteria_AnnMean,...
    NPPsum_depthCriteria_AnnMean_clim, NPPdiat_depthCriteria_AnnMean_clim, NPPdiaz_depthCriteria_AnnMean_clim, NPPsp_depthCriteria_AnnMean_clim);

load RCP8_5_2081to2100_ClimOnly.mat

outputEndCent = v2struct(yrslist, depthCriteria_clim, depthCriteria_wAnnMean, depthCriteria_wAnnMean_clim,...
    POCflux_depthCriteria_clim, POCflux_depthCriteria_AnnMean, POCflux_depthCriteria_AnnMean_clim,...
    NPPsum_depthCriteria_clim, NPPdiat_depthCriteria_clim, NPPdiaz_depthCriteria_clim, NPPsp_depthCriteria_clim,...
    NPPsum_depthCriteria_AnnMean, NPPdiat_depthCriteria_AnnMean, NPPdiaz_depthCriteria_AnnMean, NPPsp_depthCriteria_AnnMean,...
    NPPsum_depthCriteria_AnnMean_clim, NPPdiat_depthCriteria_AnnMean_clim, NPPdiaz_depthCriteria_AnnMean_clim, NPPsp_depthCriteria_AnnMean_clim);

%% Concatenate the beginning and end of century data
output = [outputToday, outputEndCent];
    
%% Plot correlation between early 21st century maximum annual MLD and change in maximum annual MLD over 21st century
% Note that this should be done with native grid, not regridded data for plotting

fignum = 1;
    figure(fignum); clf; fignum = fignum + 1;
set(gcf,'color','w')
x0=5; y0=5;
width=12; height=7;
set(gcf,'units','centimeters','position',[x0,y0,width,height])
F = 10; colormap(parula);

criteriaToCompare1 = output(1).depthCriteria_wAnnMean_clim{4}; % maximum annual MLD, 2004-2024
criteriaToCompare2 = output(2).depthCriteria_wAnnMean_clim{4}; % maximum annual MLD, 2081-2100
depthCriteriaForCorrelation = output(1).depthCriteria_wAnnMean_clim{4}; % maximum annual MLD, 2024-2024
depthCriteriaTitle = {'Maximum annual MLD, 2005-2024 (m)'};
    depthtol = 60; % don't include samples where depthCriteriaForCorrelation is less than this value

A = criteriaToCompare2 - criteriaToCompare1; A = A(:);
B = depthCriteriaForCorrelation(:);
ind = find(B > depthtol & isnan(A + B) == 0);
[rho,df,rho_sig95] = correlate(B(ind),A(ind));

dscatter(B(ind),A(ind)); hold on;
plot([1:1400],0*ones(1400,1),'k--'); hold on;
plot([1:1400],[-1:-1:-1400],'r--'); hold on;
xlabel(depthCriteriaTitle,'Fontsize',F)
ylabel({'Change in maximum annual MLD'; 'over 21st century (m)'}, 'Fontsize',F)

    figure(fignum); clf; fignum = fignum + 1;
set(gcf,'color','w'); set(gcf,'units','centimeters','position',[x0,y0,width,height])
dscatter(B(ind),A(ind)./B(ind)*100); hold on;
plot([1:1400],0*ones(1400,1),'k--'); hold on;
xlabel(depthCriteriaTitle,'Fontsize',F)
ylabel({'Percent change in maximum'; 'annual MLD over 21st century (%)'}, 'Fontsize',F)


%% Calculate fraction of NPP from diatoms
    j = 5; %Calculated at 100 m depth horizon (shouldn't matter much b/c most NPP is shallower)
fractDiatToday = output(1).NPPdiat_depthCriteria_AnnMean_clim{j}./output(1).NPPsum_depthCriteria_AnnMean_clim{j};
fractDiatEndCent = output(2).NPPdiat_depthCriteria_AnnMean_clim{j}./output(2).NPPsum_depthCriteria_AnnMean_clim{j};

    figure(fignum); clf; fignum = fignum + 1;
        cmin = -60; cmax = 60;
        C = cmocean('balance');
    subplot(311)
        imagesc(rot90(fractDiatToday)); colormap(C); colorbar;
        %caxis([cmin cmax])
        title('Fraction NPP from diatoms, 2005-2024')
    subplot(312)
        imagesc(rot90(fractDiatEndCent)); colormap(C); colorbar;
        %caxis([cmin cmax])
        title('Fraction NPP from diatoms, 2081-2100')
    subplot(313)
        imagesc(rot90(fractDiatEndCent - fractDiatToday)); colormap(C); colorbar;
        caxis([-0.3 0.3])
        title('Change in fraction NPP from diatoms, 2005-2024 to 2081-2100')

%% Calculate regional and global rates
    latint = 4;
latdiv = [-66:latint:66]; %define latitude divisions for zonal regions
numCriteria = 5; %number of depth horizons

for j = 1:2
    clear NPP_sum_glob NPP_sum_reg POCflux_sum_glob POCflux_sum_reg area_reg
    for i = 1:numCriteria
        [output(j).NPP_sum_glob(i), output(j).NPP_sum_reg(:,i), area_reg] = areaRateSum(output(j).NPPsum_depthCriteria_AnnMean_clim{i},latdiv,TLAT,TAREA,REGION_MASK);
        [output(j).POCflux_sum_glob(i), output(j).POCflux_sum_reg(:,i), area_reg] = areaRateSum(output(j).POCflux_depthCriteria_AnnMean_clim{i},latdiv,TLAT,TAREA,REGION_MASK);
    end
end

%% Plot zonal mean values and global values
cmToMeters = 1/100;

    molToPg = 12*10^-15;
    buffer = 0.1; %for plotting
    int = 3; % for plotting
    lineC = [nicecolor('mw'); nicecolor('gcck'); nicecolor('ryy'); nicecolor('k'); nicecolor('kww')];
    Lwid = 4.5-2;
    FS = 10;
for j = 1:2
    figure(fignum); clf; fignum = fignum + 1;  
set(gcf,'color','w')
x0=5;
y0=5;
width=20;
height=11.25;
set(gcf,'units','centimeters','position',[x0,y0,width,height])

for k = 1:3
if k == 1
    plotting = output(j).NPP_sum_reg./(repmat(area_reg,5,1)')/molToPg;
    titleText = {'Zonal mean NPP'};
    labelText = {'mol C m^{-2} yr^{-1}'};
    ymin = 0; ymax = 20;
    plotnum = 1;
elseif k == 2
    plotting = output(j).POCflux_sum_reg./(repmat(area_reg,5,1)')/molToPg;    %/repmat(POCflux_sum_glob,7,1); %mean annual POC flux, mol C m-2 yr-1
    titleText = {'Zonal mean POC flux'};
    labelText = {'mol C m^{-2} yr^{-1}'};
    ymin = 0; ymax = 4.8;
    plotnum = 2;
elseif k == 3
    plotting = output(j).POCflux_sum_reg./output(j).NPP_sum_reg; %mean annual e-ratio
    titleText = {'Zonal mean e-ratio'};
    labelText = {'POC flux/NPP'};
    ymin = 0; ymax = 0.4;
    plotnum = 3;
end
    subplot(3,1,plotnum)
    for i = [5,1,2,3,4]
        plot([1:length(area_reg)],plotting(:,i),'-','color',lineC(i,:),'linewidth',Lwid); hold on;
    end
    xlim([1 - buffer length(area_reg) + buffer]); xticks([2:int:length(area_reg)-1] + 0.5);
    xticklabels(latdiv([2:int:length(area_reg)-1]) + latint/2); ylim([ymin ymax]);
    ylabel(labelText,'Fontsize',FS);
    title(titleText,'Fontsize',FS);
    set(gca,'Fontsize',FS);
    if k == 3
        xlabel('Latitude','Fontsize',FS)
        legend('100 meters','Seasonal MLD','Compensation depth','Euphotic depth','Maximum annual MLD')
    end
end

end

%% Pull out seasonal cycle data from a few locations of interest

C = [nicecolor('bk'); nicecolor('ry')];
L = 3;

latstn = [60 57 40 -48]; %56, 310 - S of Greenland
lonstn = [303 336 170 160];
stnname = {'Labrador Sea','Iceland Basin','Kuroshio Extension','South Pacific'};
    nstn = length(latstn);

figure(7); clf
for k = 1:nstn %loop over all stations
% Find nearest point in CESM grid for each station of interest
    [inds(k,:),gridpt(k,:)] = findNearestPoint(TLAT,TLONG,latstn(k),lonstn(k));
mld = NaN*ones(12,2); NPP = NaN*ones(12,2); NPPdiat = NaN*ones(12,2); POCflux_100 = NaN*ones(12,2); POCflux_mldmax = NaN*ones(12,2);
for j = 1:2
    mld(:,j) = squeeze(output(j).depthCriteria_clim{1}(inds(k,1),inds(k,2),:)); %seasonal mld (12 months per grid cell)
        i = 5; %100 m depth horizon
    NPP(:,j) = squeeze(output(j).NPPsum_depthCriteria_clim{i}(inds(k,1),inds(k,2),:)); %monthly mld at given depth horizon (i)
    NPPdiat(:,j) = squeeze(output(j).NPPdiat_depthCriteria_clim{i}(inds(k,1),inds(k,2),:)); %monthly mld at given depth horizon (i)
    POCflux_100(:,j) = squeeze(output(j).POCflux_depthCriteria_clim{i}(inds(k,1),inds(k,2),:)); %monthly POC flux at given depth horizon (i)
        i = 4; %maximum annual MLD depth horizon
    POCflux_mldmax(:,j) = squeeze(output(j).POCflux_depthCriteria_clim{i}(inds(k,1),inds(k,2),:)); %monthly POC flux at given depth horizon (i)
    
    subplot(nstn,4,k)
        plot([1:24],[mld(:,j); mld(:,j)], '-', 'color', C(j,:), 'linewidth', L); hold on;
        set(gca,'ydir','reverse'); xlim([1 24]); xticks([6:6:24]); xticklabels({'Jun' 'Dec' 'Jun' 'Dec'}); title([stnname(k) 'MLD (m)']);
        ylim([0 ceil(max(max(mld))/10)*10])
    subplot(nstn,4,k + nstn)
        plot([1:24],[NPP(:,j); NPP(:,j)], '-', 'color', C(j,:), 'linewidth', L); hold on;
        xlim([1 24]); xticks([6:6:24]); xticklabels({'Jun' 'Dec' 'Jun' 'Dec'}); title('NPP (mol C m^{-2} yr^{-1})');
    subplot(nstn,4,k + nstn*3)
        plot([1:24],[NPPdiat(:,j)./NPP(:,j); NPPdiat(:,j)./NPP(:,j)], '-', 'color', C(j,:), 'linewidth', L); hold on;
        xlim([1 24]); xticks([6:6:24]); xticklabels({'Jun' 'Dec' 'Jun' 'Dec'}); title('Fraction NPP from diatoms');
    subplot(nstn,4,k + nstn*2)
        plot([1:24],[POCflux_100(:,j); POCflux_100(:,j)], '--', 'color', C(j,:), 'linewidth', L/2); hold on;
        plot([1:24],[POCflux_mldmax(:,j); POCflux_mldmax(:,j)], '-', 'color', C(j,:), 'linewidth', L); hold on;
        xlim([1 24]); xticks([6:6:24]); xticklabels({'Jun' 'Dec' 'Jun' 'Dec'}); title('POC flux (mol C m^{-2} yr^{-1})');
%     subplot(155)
%         %plot([1:24],[POCflux_100(:,j); POCflux_100(:,j)], '--', 'color', C(j,:), 'linewidth', L); hold on;
%         plot([1:24],[POCflux_mldmax(:,j); POCflux_mldmax(:,j)], '-', 'color', C(j,:), 'linewidth', L); hold on;
%         xlim([1 24]); title('POC flux at MLD_{max} (mol C m^{-2} yr^{-1})');
end

end


