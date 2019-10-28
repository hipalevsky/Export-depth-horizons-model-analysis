%% Taylor decomposition of influences on POC flux change over the 21st century
% Run after CSMBGC_RCP8_5_changeOverTime

%% Load in climatological annual mean POC flux and NPP in 3D (lat, lon, depth) for 2005-2024
load('C:\Users\palevsky\Dropbox\MATLAB\CSMBGC data processing\RCP 8_5 processed output\RCP8_5_2005to2024_ClimPOCNPPsections.mat',...
    'POCflux_all_AnnMean','NPPsum_all_AnnMean','POCflux_all_AnnMean_clim','NPPsum_all_AnnMean_clim')

%% Calculate POC flux and NPP in each year at the climatological MLDmax for 2005-2024
    TodayMLDmax_clim{1} = repmat(output(1).depthCriteria_wAnnMean_clim{4},1,1,20); %this is the climatological MLD_max for 2005-2024
POCfluxToday_TodayMLDmax_clim = depthCriteriaVals(POCflux_all_AnnMean, z_top, TodayMLDmax_clim);
NPPsumToday_TodayMLDmax_clim = depthCriteriaVals(NPPsum_all_AnnMean, z_top, TodayMLDmax_clim);

%% Validation of using climatological mean
%compare POC flux from 2005-2024 evaluated at year-specific MLDmax before
%calculating climatology with POC flux from 2005-2024 evaluated at the
%climatological MLDmax in each year
    C = cmocean('balance');      
figure(8); clf
    subplot(311)
imagesc(rot90(mean(POCfluxToday_TodayMLDmax_clim{1},3)));
colormap(C); colorbar;
caxis([0 6])
title('2005-2024 POC flux at MLD_{max}, from climatological mean MLD_{max}')
    subplot(312)
imagesc(rot90(POCflux_depthCriteria_AnnMean_clim{4}));
colormap(C); colorbar;
caxis([0 6])
title('2005-2024 POC flux at MLD_{max}, from year-specific MLD_{max}')
    subplot(313)
imagesc(rot90((POCflux_depthCriteria_AnnMean_clim{4} - ...
mean(POCfluxToday_TodayMLDmax_clim{1},3))));
colormap(C); colorbar;
caxis([-3 3])
title('Difference between POC flux at MLD_{max} from year-specific MLD_{max} vs climatological mean MLD_{max}')

%% Load in climatological annual mean POC flux and NPP in 3D (lat, lon, depth) for 2005-2024
load('C:\Users\palevsky\Dropbox\MATLAB\CSMBGC data processing\RCP 8_5 processed output\RCP8_5_2081to2100_ClimPOCNPPsections.mat',...
    'POCflux_all_AnnMean','NPPsum_all_AnnMean','POCflux_all_AnnMean_clim','NPPsum_all_AnnMean_clim')

%% Calculate mean POC depth profile for Taylor decomposition
%Reshape as list of depth profiles in all grid cells
[len, wid, h] = size(POCflux_all_AnnMean_clim);
POCflux_reshaped = NaN*ones(len*wid, h);
area_weighting = NaN*ones(len*wid, 1);
counter = 1;
for i = 1:len
    for j = 1:wid
        POCflux_reshaped(counter, :) = squeeze(POCflux_all_AnnMean_clim(i,j,:));
        area_weighting(counter, :) = TAREA(i,j);
        counter = counter + 1;
    end
end

POCflux_globmean = nanmean(POCflux_reshaped.*repmat(area_weighting, 1, 60))/nanmean(area_weighting);
POCflux_globmean_noweight = nanmean(POCflux_reshaped);
POCflux_globmean_interp = interp1(z_top, POCflux_globmean, [1:1400]);

%% Plot global mean POC flux profile
figure(10); clf
    set(gcf,'color','w')
    x0=5; y0=5;
    width=7; height=7;
    set(gcf,'units','centimeters','position',[x0,y0,width,height])
    L = 2; M = 8;
    
plot(POCflux_globmean_interp, [1:1400], '-', 'color', nicecolor('bcw'), 'linewidth', L); hold on;
plot(POCflux_globmean, z_top, 'k.','markersize',M); hold on;
plot([0 2.2], [100 100], 'r--', 'linewidth', L/1.5); hold on;
set(gca,'YDir','reverse'); 
ylim([0 1000]); xlim([0 2.2])
title('Global mean POC flux profile')
ylabel('Depth (m)')
xlabel('mol C m^{-2} yr^{-1}')

%% Calculate first-order Taylor decomposition of POC flux changes (see Laufkotter et al. 2016, Eqn 18, Fig 5)
% Note that these initial calculations start with the climatological annual
% mean, but could go back to calculate for each individual year before
% taking mean

%Calculate Taylor residual for dFlux/dz * dMLD/dt based on the global mean
%flux profile from the beginning of the century
MLDmax_endcent = ceil(output(2).depthCriteria_wAnnMean_clim{4}); %take ceiling to get an ID to use with POCflux_globmean_interp
MLDmax_begcent = ceil(output(1).depthCriteria_wAnnMean_clim{4}); %take ceiling to get an ID to use with POCflux_globmean_interp

Taylor_MLDchange = NaN*ones(len,wid); %initialize with NaNs
for i = 1:len
    for j = 1:wid
        if isnan(MLDmax_endcent(i,j) + MLDmax_begcent(i,j)) == 0
            Taylor_MLDchange(i,j) = POCflux_globmean_interp(MLDmax_endcent(i,j)) - POCflux_globmean_interp(MLDmax_begcent(i,j));
        end
    end
end

j = 4; %Calculate change in EP at MLD_max
dEPdt_MLDmax = output(2).POCflux_depthCriteria_AnnMean_clim{j} - output(1).POCflux_depthCriteria_AnnMean_clim{j};
deRatiodt_MLDmax = output(2).POCflux_depthCriteria_AnnMean_clim{j}./output(2).NPPsum_depthCriteria_AnnMean_clim{j} -...
        output(1).POCflux_depthCriteria_AnnMean_clim{j}./output(1).NPPsum_depthCriteria_AnnMean_clim{j};

j = 5; %100 m - calculate baseline values at 100 m and then add in dMLD/dt term from above
    dEPdt = output(2).POCflux_depthCriteria_AnnMean_clim{j} - output(1).POCflux_depthCriteria_AnnMean_clim{j};
    dNPPdt = output(2).NPPsum_depthCriteria_AnnMean_clim{j} - output(1).NPPsum_depthCriteria_AnnMean_clim{j};
    deRatiodt = output(2).POCflux_depthCriteria_AnnMean_clim{j}./output(2).NPPsum_depthCriteria_AnnMean_clim{j} -...
        output(1).POCflux_depthCriteria_AnnMean_clim{j}./output(1).NPPsum_depthCriteria_AnnMean_clim{j};

%Calculate Taylor terms using early century clim values for non-partial
%derivative terms in equation
    Taylor_NPPTerm = dNPPdt.*(output(1).POCflux_depthCriteria_AnnMean_clim{j}./output(1).NPPsum_depthCriteria_AnnMean_clim{j});
    Taylor_eRatioTerm = deRatiodt.*(output(1).NPPsum_depthCriteria_AnnMean_clim{j});
    Taylor_eRatioTerm_MLDmax = deRatiodt_MLDmax.*(output(1).NPPsum_depthCriteria_AnnMean_clim{4});
    Taylor_Residual_100m = dEPdt - Taylor_NPPTerm - Taylor_eRatioTerm;
    Taylor_Residual_MLDmax = dEPdt_MLDmax - Taylor_NPPTerm - Taylor_eRatioTerm_MLDmax; %Scott suggested adding Taylor_MLDchange, but residual still doesn't work out right
    Taylor_Residual_MLDmax_wMLDchange = dEPdt_MLDmax - Taylor_NPPTerm - Taylor_eRatioTerm - Taylor_MLDchange; 
    
    figure(fignum); clf; fignum = fignum + 1;
        cmin = -100; cmax = 100;
        C = cmocean('balance');
        j = 4; %For MLD_max
    subplot(511)
        imagesc(rot90(100*dEPdt_MLDmax./output(1).POCflux_depthCriteria_AnnMean_clim{j})); colormap(C); colorbar;
        caxis([cmin cmax])
        title('Percent change in POC flux, from 2005-2024 to 2081-2100')
    subplot(512)
        imagesc(rot90(100*Taylor_NPPTerm./output(1).POCflux_depthCriteria_AnnMean_clim{j})); colormap(C); colorbar;
        caxis([cmin cmax])
        title('Percent change due to NPP (\deltaNPP/\deltat x e-ratio)')
    subplot(513)
        imagesc(rot90(100*Taylor_eRatioTerm_MLDmax./output(1).POCflux_depthCriteria_AnnMean_clim{j})); colormap(C); colorbar;
        caxis([cmin cmax])
        title('Percent change due to e-ratio (\deltae-ratio/\deltat x NPP)')
    subplot(514)
        imagesc(rot90(100*Taylor_MLDchange./output(1).POCflux_depthCriteria_AnnMean_clim{j})); colormap(C); colorbar;
        caxis([cmin cmax])
        title('Percent change due to MLD_{max} (\deltaPOC flux/\deltaz x \deltaMLD_{max})')
    subplot(515)
        imagesc(rot90(100*Taylor_Residual_MLDmax./output(1).POCflux_depthCriteria_AnnMean_clim{j})); colormap(C); colorbar;
        caxis([cmin cmax])
        title('Residual of 1st order Taylor decomposition above (%)')
    