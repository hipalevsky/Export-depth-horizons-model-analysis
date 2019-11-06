%% Taylor decomposition of influences on POC flux change over the 21st century
% Run after CSMBGC_RCP8_5_changeOverTime

% This script calculates the Taylor decomposition at both the fixed 100 m depth
% horizon (as was done in Laufkotter et al 2016) and at the MLD_max
% horizon, following three different methods:
% 1) Wrapping all of the MLD_max change into the E-Ratio term,
% 2) Calculating NPP and E-Ratio terms at the fixed 100 m depth horizon and 
    % then adding a third term to account for POC flux change due to change
    % in MLD, using the global mean POC flux curve to estimate the influence,
% 3) Calculating NPP and E-Ratio terms at the climatological MLD_max from 
    % the beginning of the century (2005-2024) and adding a third terms to 
    % account for POC flux change due to change in MLD by taking the 
    % difference between the actual POC flux at the 2081-2100 MLD_max and
    % the POC flux from those same years but instead determined at the
    % climatological MLD_max from 2005-2024.

%% Load in climatological annual mean POC flux and NPP in 3D (lat, lon, depth) for 2005-2024
load('C:\Users\palevsky\Dropbox\MATLAB\CSMBGC data processing\RCP 8_5 processed output\RCP8_5_2005to2024_ClimPOCNPPsections.mat',...
    'POCflux_all_AnnMean','NPPsum_all_AnnMean','POCflux_all_AnnMean_clim','NPPsum_all_AnnMean_clim')

%% Calculate POC flux and NPP in each year
%at the climatological MLDmax for 2005-2024 and at the climatological MLDmax for 2081-2100
    MLDmax_clim{1} = repmat(output(1).depthCriteria_wAnnMean_clim{4},1,1,20); %this is the climatological MLD_max for 2005-2024
    MLDmax_clim{2} = repmat(output(2).depthCriteria_wAnnMean_clim{4},1,1,20); %this is the climatological MLD_max for 2081-2100

%Evaluate today's fluxes at the climatological MLDmax for both 2005-2024 and 2018-2100
POCfluxToday_MLDmax_clim = depthCriteriaVals(POCflux_all_AnnMean, z_top, MLDmax_clim);
NPPsumToday_MLDmax_clim = depthCriteriaVals(NPPsum_all_AnnMean, z_top, MLDmax_clim);

%% Validation of using climatological mean
%compare POC flux from 2005-2024 evaluated at year-specific MLDmax before
%calculating climatology with POC flux from 2005-2024 evaluated at the
%climatological MLDmax in each year
    C = cmocean('balance');      
figure(8); clf
    subplot(311)
imagesc(rot90(mean(POCfluxToday_MLDmax_clim{1},3)));
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
mean(POCfluxToday_MLDmax_clim{1},3))));
colormap(C); colorbar;
caxis([-3 3])
title('Difference between POC flux at MLD_{max} from year-specific MLD_{max} vs climatological mean MLD_{max}')

%% Load in climatological annual mean POC flux and NPP in 3D (lat, lon, depth) for 2005-2024
load('C:\Users\palevsky\Dropbox\MATLAB\CSMBGC data processing\RCP 8_5 processed output\RCP8_5_2081to2100_ClimPOCNPPsections.mat',...
    'POCflux_all_AnnMean','NPPsum_all_AnnMean','POCflux_all_AnnMean_clim','NPPsum_all_AnnMean_clim')

%% Calculate POC flux and NPP in each year from 2081-2100 at the climatological MLDmax for 2005-2024 and for 2081-2100
POCfluxEndCent_MLDmax_clim = depthCriteriaVals(POCflux_all_AnnMean, z_top, MLDmax_clim);
NPPsumEndCent_MLDmax_clim = depthCriteriaVals(NPPsum_all_AnnMean, z_top, MLDmax_clim);

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

%Calculate Taylor decomposition for MLD_max using climatological MLD_max
%values from beginning and end of century as extra depth criteria
dNPPdt_MLDmaxBegCent = mean(NPPsumEndCent_MLDmax_clim{1},3) - mean(NPPsumToday_MLDmax_clim{1},3);
eRatio_MLDmaxBegCent = mean((POCfluxToday_MLDmax_clim{1}./NPPsumToday_MLDmax_clim{1}),3);
Taylor_NPPTerm_MLDmaxBegCent = dNPPdt_MLDmaxBegCent.*eRatio_MLDmaxBegCent;

dERatiodt_MLDmaxBegCent = mean((POCfluxEndCent_MLDmax_clim{1}./NPPsumEndCent_MLDmax_clim{1}),3) - ...
    mean((POCfluxToday_MLDmax_clim{1}./NPPsumToday_MLDmax_clim{1}),3);
NPP_MLDmaxBegCent = mean(NPPsumToday_MLDmax_clim{1},3);
Taylor_eRatioTerm_MLDmaxBegCent = dERatiodt_MLDmaxBegCent.*NPP_MLDmaxBegCent;

%Calculated with climatological mean MLD_max depth horizons for both
%beginning and end of century
    %Taylor_dEPdt_MLDmaxChange = mean(POCfluxEndCent_MLDmax_clim{2},3) - mean(POCfluxEndCent_MLDmax_clim{1},3); 
%Calculated with climatological MLD_max depth horizon from beginning of
%century, but depth-specific values for end of century - BETTER METHOD
Taylor_dEPdt_MLDmaxChange = output(2).POCflux_depthCriteria_AnnMean_clim{4} - mean(POCfluxEndCent_MLDmax_clim{1},3);

Taylor_Residual_MLDmax_BegVsEndCentCriteria = dEPdt_MLDmax - Taylor_NPPTerm_MLDmaxBegCent -...
    Taylor_eRatioTerm_MLDmaxBegCent - Taylor_dEPdt_MLDmaxChange;
    
%% Visualize Taylor Decomposition terms for the two different ways of adding a 3rd term in the MLD_max decomposition

    figure(fignum); clf; fignum = fignum + 1;
        cmin = -1.5; cmax = 1.5;
        C = cmocean('balance');
        j = 4; %For MLD_max
    subplot(511)
        imagesc(rot90(dEPdt_MLDmax)); colormap(C); colorbar;
        caxis([cmin cmax])
        title('Change in POC flux at MLD_{max}, from 2005-2024 to 2081-2100')
    subplot(512)
        imagesc(rot90(Taylor_NPPTerm_MLDmaxBegCent)); colormap(C); colorbar;
        caxis([cmin cmax])
        title('Change due to NPP at MLD_{max, 2005-2024} (\deltaNPP/\deltat x e-ratio)')
    subplot(513)
        imagesc(rot90(Taylor_eRatioTerm_MLDmaxBegCent)); colormap(C); colorbar;
        caxis([cmin cmax])
        title('Change due to e-ratio at MLD_{max, 2005-2024} (\deltae-ratio/\deltat x NPP)')
    subplot(514)
        imagesc(rot90(Taylor_dEPdt_MLDmaxChange)); colormap(C); colorbar;
        caxis([cmin cmax])
        title('Change due to change in MLD_{max} depth criterion (\deltaPOC flux/\deltaz x \deltaMLD_{max})')
    subplot(515)
        imagesc(rot90(Taylor_Residual_MLDmax_BegVsEndCentCriteria)); colormap(C); colorbar;
        caxis([cmin cmax])
        title('Residual of 1st order Taylor decomposition above')

% Visualize 2nd method of taylor decomposition w/ beginning and end of 21st century-specific MLD_max horizons        
    figure(fignum); clf; fignum = fignum + 1;
        cmin = -1.5; cmax = 1.5;
        C = cmocean('balance');
        j = 4; %For MLD_max
    subplot(511)
        imagesc(rot90(dEPdt_MLDmax)); colormap(C); colorbar;
        caxis([cmin cmax])
        title('Change in POC flux, from 2005-2024 to 2081-2100')
    subplot(512)
        imagesc(rot90(Taylor_NPPTerm)); colormap(C); colorbar;
        caxis([cmin cmax])
        title('Change due to NPP (\deltaNPP/\deltat x e-ratio)')
    subplot(513)
        imagesc(rot90(Taylor_eRatioTerm_MLDmax)); colormap(C); colorbar;
        caxis([cmin cmax])
        title('Percent change due to e-ratio (\deltae-ratio/\deltat x NPP)')
    subplot(514)
        imagesc(rot90(Taylor_MLDchange)); colormap(C); colorbar;
        caxis([cmin cmax])
        title('Change due to MLD_{max} (\deltaPOC flux/\deltaz x \deltaMLD_{max})')
    subplot(515)
        imagesc(rot90(Taylor_Residual_MLDmax)); colormap(C); colorbar;
        caxis([cmin cmax])
        title('Residual of 1st order Taylor decomposition above')

% Compare the residuals with the total change
   figure(fignum); clf; fignum = fignum + 1;
   cmin = -1.5; cmax = 1.5;
       subplot(311)
        imagesc(rot90(dEPdt_MLDmax)); colormap(C); colorbar;
        caxis([cmin cmax])
        title('Total dEPdt at MLD_max')
        
      subplot(312)
        imagesc(rot90(Taylor_Residual_MLDmax_BegVsEndCentCriteria)); colormap(C); colorbar;
        caxis([cmin cmax])
        title('Residual of 1st order Taylor decomposition - new method')
        
      subplot(313)
        imagesc(rot90(Taylor_Residual_MLDmax)); colormap(C); colorbar;
        caxis([cmin cmax])
        title('Residual of 1st order Taylor decomposition - old method')
        
%% Clear variables to free up memory
clear POCflux_all_AnnMean NPPsum_all_AnnMean POCflux_all_AnnMean_clim NPPsum_all_AnnMean_clim
    