load RCP8_5_GriddedOutput.mat %saved in CSMBGC_RCP8_5_changeOverTime_regrid
%note that workflow to this point is:
% 1) CSMBGC_RCP8_5_changeOverTime
% 2) CSMBGC_RCP8_5_TaylorDecomposition
% 3) CSMBGC_RCP8_5_changeOverTime_regrid

%% Time series site locations to plot on maps
latstn = [60 57 40 -48]; %56, 310 - S of Greenland
lonstn = [303 336 170 160];
stnname = {'Labrador Sea','Iceland Basin','Kuroshio Extension','South Pacific'};

%% Compare MLD_max for beginning and end of century

F = 10; %fontsize

%Maximum annual MLD at beginning and end of century
    figure(1); clf;
set(gcf,'color','w')
x0=2;
y0=2;
width=16;
height=14;
set(gcf,'units','centimeters','position',[x0,y0,width,height])

latminplot = -75; latmaxplot = 75;
lonminplot = 10; lonmaxplot = 390;
cmin = 0; cmax = 380; cint = 20;
C = cmocean('deep');
set(0,'defaultAxesFontSize',F)

for i = 1:2
    subplot(2,1,i)
    m_proj('miller','lat',[latminplot latmaxplot],'lon',[lonminplot lonmaxplot])
    m_contourf([glon glon + 360], glat, repmat(outputGrid(i).depthCriteria_wAnnMean_clim_grid(:,:,4),1,2),[cmin: cint: cmax],'linecolor','none'); hold on;
    colormap(C); caxis([cmin cmax]); hc = colorbar('eastoutside');
    m_grid('box','fancy'); m_coast('patch',nicecolor('wwk'));
    if i == 1
        title('a) Climatological MLD_{max}, 2005-2024 (m)')
    elseif i == 2
        title('b) Climatological MLD_{max}, 2081-2100 (m)')
    end
end
%%
%Difference between maximum annual MLD at beginning and end of century
    figure(2); clf;
set(gcf,'color','w')
x0=5;
y0=5;
width=16;
height=7;
set(gcf,'units','centimeters','position',[x0,y0,width,height])

cmin = -140; cmax = 140; cint = 10;
C = cmocean('balance');

%Create a plotting variable with changes out of min bounds filled with min
%of colorbar (for nicer plot)
MLDmax_change_21stCent_plot = outputGrid(2).depthCriteria_wAnnMean_clim_grid(:,:,4) - ...
        outputGrid(1).depthCriteria_wAnnMean_clim_grid(:,:,4);
for i = 1:length(glat)
    indminout = find(MLDmax_change_21stCent_plot(i,:) < cmin);
    MLDmax_change_21stCent_plot(i,indminout) = cmin;
end

    m_proj('miller','lat',[latminplot latmaxplot],'lon',[lonminplot lonmaxplot])
    m_contourf([glon glon + 360], glat, repmat(MLDmax_change_21stCent_plot,1,2),[cmin: cint: cmax],'linecolor','none'); hold on;
    m_plot(lonstn(1:2), latstn(1:2), 'ko','markersize',8,'markerfacecolor','y'); hold on;
    colormap(C); caxis([cmin cmax]); hc = colorbar('eastoutside');
    m_grid('box','fancy'); m_coast('patch',nicecolor('wwk'));
    title('c) Change in MLD_{max} over the 21st century (m)')
    
%% Percent Difference between maximum annual MLD at beginning and end of century
% Just to show that places with deep initial MLDs and big absolute changes
% are also places with large % changes
    figure(3); clf;
set(gcf,'color','w')
x0=5;
y0=5;
width=16;
height=7;
set(gcf,'units','centimeters','position',[x0,y0,width,height])

cmin = -100; cmax = 100; cint = 10;
C = cmocean('balance');

%Create a plotting variable with changes out of min bounds filled with min
%of colorbar (for nicer plot)
MLDmax_percentchange_21stCent_plot = 100*(outputGrid(2).depthCriteria_wAnnMean_clim_grid(:,:,4) - ...
        outputGrid(1).depthCriteria_wAnnMean_clim_grid(:,:,4))./outputGrid(1).depthCriteria_wAnnMean_clim_grid(:,:,4);

    m_proj('miller','lat',[latminplot latmaxplot],'lon',[lonminplot lonmaxplot])
    m_contourf([glon glon + 360], glat, repmat(MLDmax_percentchange_21stCent_plot,1,2),[cmin: cint: cmax],'linecolor','none'); hold on;
    m_plot(lonstn, latstn, 'ko','markersize',6,'markerfacecolor','y'); hold on;
    colormap(C); caxis([cmin cmax]); hc = colorbar('eastoutside');
    m_grid('box','fancy'); m_coast('patch',nicecolor('wwk'));
    title('Percent change in climatological MLD_{max} over 21st century (%)')
    
%% Compare POC flux change over 21st century at 100 m and at MLD_max

    figure(4); clf;
set(gcf,'color','w')
x0=0;
y0=0;
width=16;
height=21;
set(gcf,'units','centimeters','position',[x0,y0,width,height])

latminplot = -75; latmaxplot = 75;
lonminplot = 10; lonmaxplot = 390;
cmin = -2; cmax = 2; cint = 0.2;
C = cmocean('balance');
set(0,'defaultAxesFontSize',F)

for i = 1:2
    subplot(3,1,i)
    m_proj('miller','lat',[latminplot latmaxplot],'lon',[lonminplot lonmaxplot])
    m_contourf([glon glon + 360], glat, repmat((outputGrid(2).POCflux_depthCriteria_AnnMean_clim_grid(:,:,6 - i) - ...
        outputGrid(1).POCflux_depthCriteria_AnnMean_clim_grid(:,:,6 - i)),1,2),[cmin: cint: cmax],'linecolor','none'); hold on;
    colormap(C); caxis([cmin cmax]); hc = colorbar('eastoutside');
    m_grid('box','fancy'); m_coast('patch',nicecolor('wwk'));
    if i == 1
        title('a) POC flux change evaluated at 100 m')
    elseif i == 2
        title('b) POC flux change evaluated at MLD_{max}')
    end
end

subplot(3,1,3)
%Create a plotting variable for difference in change between mldMax and 100 m depth horizons
POCflux_changeDifference_21stCent_plot = -(outputGrid(2).POCflux_depthCriteria_AnnMean_clim_grid(:,:,5) - ...
        outputGrid(1).POCflux_depthCriteria_AnnMean_clim_grid(:,:,5)) + (outputGrid(2).POCflux_depthCriteria_AnnMean_clim_grid(:,:,4) - ...
        outputGrid(1).POCflux_depthCriteria_AnnMean_clim_grid(:,:,4));
cmin = -1; cmax = 1; cint = 0.1;

    m_proj('miller','lat',[latminplot latmaxplot],'lon',[lonminplot lonmaxplot])
    m_contourf([glon glon + 360], glat, repmat(POCflux_changeDifference_21stCent_plot,1,2),[cmin: cint: cmax],'linecolor','none'); hold on;
    colormap(C); caxis([cmin cmax]); hc = colorbar('eastoutside');
    m_grid('box','fancy'); m_coast('patch',nicecolor('wwk'));
    title('c) Difference in POC flux change, MLD_{max} - 100 m')
    
%% Compare e-ratio change over 21st century at 100 m and at MLD_max

    figure(5); clf;
set(gcf,'color','w')
x0=0;
y0=0;
width=16;
height=21;
set(gcf,'units','centimeters','position',[x0,y0,width,height])

latminplot = -75; latmaxplot = 75;
lonminplot = 10; lonmaxplot = 390;
cmin = -0.14; cmax = 0.14; cint = 0.01;
C = cmocean('balance');
set(0,'defaultAxesFontSize',F)

for i = 1:2
    subplot(3,1,i)
    m_proj('miller','lat',[latminplot latmaxplot],'lon',[lonminplot lonmaxplot])
    m_contourf([glon glon + 360], glat, repmat((outputGrid(2).POCflux_depthCriteria_AnnMean_clim_grid(:,:,6 - i)./outputGrid(2).NPPsum_depthCriteria_AnnMean_clim_grid(:,:,6 - i) - ...
        outputGrid(1).POCflux_depthCriteria_AnnMean_clim_grid(:,:,6 - i)./outputGrid(1).NPPsum_depthCriteria_AnnMean_clim_grid(:,:,6 - i)),1,2),...
        [cmin: cint: cmax],'linecolor','none'); hold on;
    colormap(C); caxis([cmin cmax]); hc = colorbar('eastoutside');
    m_grid('box','fancy'); m_coast('patch',nicecolor('wwk'));
    if i == 1
        title('a) e-ratio change evaluated at 100 m')
    elseif i == 2
        title('b) e-ratio change evaluated at MLD_{max}')
    end
end

subplot(3,1,3)
%Create a plotting variable for difference in change between mldMax and 100 m depth horizons
eratio_changeDifference_21stCent_plot = -(outputGrid(2).POCflux_depthCriteria_AnnMean_clim_grid(:,:,5)./outputGrid(2).NPPsum_depthCriteria_AnnMean_clim_grid(:,:,5) - ...
        outputGrid(1).POCflux_depthCriteria_AnnMean_clim_grid(:,:,5)./outputGrid(1).NPPsum_depthCriteria_AnnMean_clim_grid(:,:,5)) + ...
        (outputGrid(2).POCflux_depthCriteria_AnnMean_clim_grid(:,:,4)./outputGrid(2).NPPsum_depthCriteria_AnnMean_clim_grid(:,:,4) - ...
        outputGrid(1).POCflux_depthCriteria_AnnMean_clim_grid(:,:,4)./outputGrid(1).NPPsum_depthCriteria_AnnMean_clim_grid(:,:,4));
cmin = -0.58; cmax = 0.58; cint = 0.05;

    m_proj('miller','lat',[latminplot latmaxplot],'lon',[lonminplot lonmaxplot])
    m_contourf([glon glon + 360], glat, repmat(POCflux_changeDifference_21stCent_plot,1,2),[cmin: cint: cmax],'linecolor','none'); hold on;
    colormap(C); caxis([cmin cmax]); hc = colorbar('eastoutside');
    m_grid('box','fancy'); m_coast('patch',nicecolor('wwk'));
    title('c) Difference in e-ratio change, MLD_{max} - 100 m')
    
%% Calculate diatom fraction of NPP at beginning and of century, and change over time

    figure(6); clf;
set(gcf,'color','w')
x0=0;
y0=0;
width=16;
height=14;
set(gcf,'units','centimeters','position',[x0,y0,width,height])

C = cmocean('deep');
set(0,'defaultAxesFontSize',F)
    j = 5; %set to use 100 m depth horizon
fractDiatToday = outputGrid(1).NPPdiat_depthCriteria_AnnMean_clim_grid(:,:,j)./outputGrid(1).NPPsum_depthCriteria_AnnMean_clim_grid(:,:,j);
fractDiatEndCent = outputGrid(2).NPPdiat_depthCriteria_AnnMean_clim_grid(:,:,j)./outputGrid(2).NPPsum_depthCriteria_AnnMean_clim_grid(:,:,j);
cmin = 0; cmax = 100; cint = 1;

subplot(2,1,1)
    m_proj('miller','lat',[latminplot latmaxplot],'lon',[lonminplot lonmaxplot])
    m_contourf([glon glon + 360], glat, repmat(100*fractDiatToday,1,2), [cmin: cint: cmax],'linecolor','none'); hold on;
    colormap(C); caxis([cmin cmax]); hc = colorbar('eastoutside');
    m_grid('box','fancy'); m_coast('patch',nicecolor('wwk'));
    title('Percent NPP from diatoms, 2005-2024 (%)')
    
subplot(2,1,2)
    m_proj('miller','lat',[latminplot latmaxplot],'lon',[lonminplot lonmaxplot])
    m_contourf([glon glon + 360], glat, repmat(100*fractDiatEndCent,1,2), [cmin: cint: cmax],'linecolor','none'); hold on;
    colormap(C); caxis([cmin cmax]); hc = colorbar('eastoutside');
    m_grid('box','fancy'); m_coast('patch',nicecolor('wwk'));
    title('Percent NPP from diatoms, 2081-2100 (%)')

   
    figure(7); clf;
set(gcf,'color','w')
x0=0;
y0=0;
width=16;
height=7;
set(gcf,'units','centimeters','position',[x0,y0,width,height])
C = cmocean('balance');

    cmin = -38; cmax = 38; cint = 1;
fractDiatChange = (fractDiatEndCent - fractDiatToday)*100;
for i = 1:length(glat)
    indminout = find(fractDiatChange(i,:) < cmin);
    fractDiatChange(i,indminout) = cmin;
end
    m_proj('miller','lat',[latminplot latmaxplot],'lon',[lonminplot lonmaxplot])
    m_contourf([glon glon + 360], glat, repmat(fractDiatChange,1,2), [cmin: cint: cmax],'linecolor','none'); hold on;
    colormap(C); caxis([cmin cmax]); hc = colorbar('eastoutside');
    m_grid('box','fancy'); m_coast('patch',nicecolor('wwk'));
    title('Change in percent NPP from diatoms, 2005-2024 to 2081-2100')

%% Taylor decomposition of change in POC flux over time

    %Decomposition at MLDmax, without extra term
figure(8); clf;
set(gcf,'color','w')
x0=0;
y0=0;
width=34;
height=14;
set(gcf,'units','centimeters','position',[x0,y0,width,height])
    cmin = -75; cmax = 75;
    C = cmocean('balance');
    
    NormalizeData = outputGrid(1).POCflux_depthCriteria_AnnMean_clim_grid(:,:,4); %Annual mean POC flux from beginning of century

subplot(221)
    m_proj('miller','lat',[latminplot latmaxplot],'lon',[lonminplot lonmaxplot])
    m_contourf([glon glon + 360], glat, repmat(100*(TaylorGrid.dEPdt_MLDmax./NormalizeData),1,2), [cmin: cint: cmax],'linecolor','none'); hold on;
    colormap(C); caxis([cmin cmax]); hc = colorbar('eastoutside');
    m_grid('box','fancy'); m_coast('patch',nicecolor('wwk'));
    title('Percent change in POC flux at MLD_{max}, from 2005-2024 to 2081-2100')
subplot(222)
    m_proj('miller','lat',[latminplot latmaxplot],'lon',[lonminplot lonmaxplot])
    m_contourf([glon glon + 360], glat, repmat(100*(TaylorGrid.NPPterm./NormalizeData),1,2), [cmin: cint: cmax],'linecolor','none'); hold on;
    colormap(C); caxis([cmin cmax]); hc = colorbar('eastoutside');
    m_grid('box','fancy'); m_coast('patch',nicecolor('wwk'));
    title('Percent change due to NPP (\deltaNPP/\deltat x e-ratio)')
subplot(223)
    m_proj('miller','lat',[latminplot latmaxplot],'lon',[lonminplot lonmaxplot])
    m_contourf([glon glon + 360], glat, repmat(100*(TaylorGrid.eRatioTerm_MLDmax./NormalizeData),1,2), [cmin: cint: cmax],'linecolor','none'); hold on;
    colormap(C); caxis([cmin cmax]); hc = colorbar('eastoutside');
    m_grid('box','fancy'); m_coast('patch',nicecolor('wwk'));
    title('Percent change due to e-ratio (\deltae-ratio/\deltat x NPP)')
subplot(224)
    m_proj('miller','lat',[latminplot latmaxplot],'lon',[lonminplot lonmaxplot])
    m_contourf([glon glon + 360], glat, repmat(100*(TaylorGrid.Residual_MLDmax./NormalizeData),1,2), [cmin: cint: cmax],'linecolor','none'); hold on;
    colormap(C); caxis([cmin cmax]); hc = colorbar('eastoutside');
    m_grid('box','fancy'); m_coast('patch',nicecolor('wwk'));
    title('Residual of Taylor decomposition for change in MLD_{max} POC flux')

  %%  Taylor decomposition at MLDmax with extra term calculated using beginning of century MLDmax depth criterion
figure(9); clf;
set(gcf,'color','w')
x0=0;
y0=0;
width=34;
height=22;
set(gcf,'units','centimeters','position',[x0,y0,width,height])
    %cmin = -75; cmax = 75;
    cmin = -2; cmax = 2; cint = 0.02;
    C = cmocean('balance');
    
    NormalizeData = 100; %outputGrid(1).POCflux_depthCriteria_AnnMean_clim_grid(:,:,4); %Annual mean POC flux from beginning of century

subplot(3,2,[1:2])
    m_proj('miller','lat',[latminplot latmaxplot],'lon',[lonminplot lonmaxplot])
    m_contourf([glon glon + 360], glat, repmat(100*(TaylorGrid.dEPdt_MLDmax./NormalizeData),1,2), [cmin: cint: cmax],'linecolor','none'); hold on;
    colormap(C); caxis([cmin cmax]); hc = colorbar('eastoutside'); ylabel(hc,'mol C m^{-2} yr^{-1}')
    m_grid('box','fancy'); m_coast('patch',nicecolor('wwk'));
    title('a) Total POC flux change evaluated at MLD_{max}')
subplot(323)
    m_proj('miller','lat',[latminplot latmaxplot],'lon',[lonminplot lonmaxplot])
    m_contourf([glon glon + 360], glat, repmat(100*(TaylorGrid.NPPTerm_MLDmaxBegCent./NormalizeData),1,2), [cmin: cint: cmax],'linecolor','none'); hold on;
    colormap(C); caxis([cmin cmax]); hc = colorbar('eastoutside'); ylabel(hc,'mol C m^{-2} yr^{-1}')
    m_grid('box','fancy'); m_coast('patch',nicecolor('wwk'));
    title('b) NPP term')
subplot(324)
    m_proj('miller','lat',[latminplot latmaxplot],'lon',[lonminplot lonmaxplot])
    m_contourf([glon glon + 360], glat, repmat(100*(TaylorGrid.eRatioTerm_MLDmaxBegCent./NormalizeData),1,2), [cmin: cint: cmax],'linecolor','none'); hold on;
    colormap(C); caxis([cmin cmax]); hc = colorbar('eastoutside'); ylabel(hc,'mol C m^{-2} yr^{-1}')
    m_grid('box','fancy'); m_coast('patch',nicecolor('wwk'));
    title('c) e-ratio term')
subplot(325)
    m_proj('miller','lat',[latminplot latmaxplot],'lon',[lonminplot lonmaxplot])
    m_contourf([glon glon + 360], glat, repmat(100*(TaylorGrid.dEPdt_MLDmaxChange./NormalizeData),1,2), [cmin: cint: cmax],'linecolor','none'); hold on;
    colormap(C); caxis([cmin cmax]); hc = colorbar('eastoutside'); ylabel(hc,'mol C m^{-2} yr^{-1}')
    m_grid('box','fancy'); m_coast('patch',nicecolor('wwk'));
    title('d) \DeltaMLD_{max} term')
subplot(326)
    m_proj('miller','lat',[latminplot latmaxplot],'lon',[lonminplot lonmaxplot])
    m_contourf([glon glon + 360], glat, repmat(100*(TaylorGrid.Residual_MLDmax_BegVsEndCentCriteria./NormalizeData),1,2), [cmin: cint: cmax],'linecolor','none'); hold on;
    colormap(C); caxis([cmin cmax]); hc = colorbar('eastoutside'); ylabel(hc,'mol C m^{-2} yr^{-1}')
    m_grid('box','fancy'); m_coast('patch',nicecolor('wwk'));
    title('e) Residual')
    
  %%  Taylor decomposition at MLDmax with extra term calculated using POC flux profile
figure(90); clf;
set(gcf,'color','w')
x0=0;
y0=0;
width=34;
height=22;
set(gcf,'units','centimeters','position',[x0,y0,width,height])
    cmin = -75; cmax = 75;
    C = cmocean('balance');
    
    NormalizeData = outputGrid(1).POCflux_depthCriteria_AnnMean_clim_grid(:,:,4); %Annual mean POC flux from beginning of century

subplot(321)
    m_proj('miller','lat',[latminplot latmaxplot],'lon',[lonminplot lonmaxplot])
    m_contourf([glon glon + 360], glat, repmat(100*(TaylorGrid.dEPdt_MLDmax./NormalizeData),1,2), [cmin: cint: cmax],'linecolor','none'); hold on;
    colormap(C); caxis([cmin cmax]); hc = colorbar('eastoutside'); 
    m_grid('box','fancy'); m_coast('patch',nicecolor('wwk'));
    title('Percent change in POC flux at MLD_{max}, from 2005-2024 to 2081-2100')
subplot(322)
    m_proj('miller','lat',[latminplot latmaxplot],'lon',[lonminplot lonmaxplot])
    m_contourf([glon glon + 360], glat, repmat(100*(TaylorGrid.NPPterm./NormalizeData),1,2), [cmin: cint: cmax],'linecolor','none'); hold on;
    colormap(C); caxis([cmin cmax]); hc = colorbar('eastoutside');
    m_grid('box','fancy'); m_coast('patch',nicecolor('wwk'));
    title('Percent change due to NPP (\deltaNPP/\deltat x e-ratio)')
subplot(323)
    m_proj('miller','lat',[latminplot latmaxplot],'lon',[lonminplot lonmaxplot])
    m_contourf([glon glon + 360], glat, repmat(100*(TaylorGrid.eRatioTerm_100m./NormalizeData),1,2), [cmin: cint: cmax],'linecolor','none'); hold on;
    colormap(C); caxis([cmin cmax]); hc = colorbar('eastoutside');
    m_grid('box','fancy'); m_coast('patch',nicecolor('wwk'));
    title('Percent change due to e-ratio at 100m (\deltae-ratio_{100m}/\deltat x NPP)')
subplot(324)
    m_proj('miller','lat',[latminplot latmaxplot],'lon',[lonminplot lonmaxplot])
    m_contourf([glon glon + 360], glat, repmat(100*(TaylorGrid.MLDchange./NormalizeData),1,2), [cmin: cint: cmax],'linecolor','none'); hold on;
    colormap(C); caxis([cmin cmax]); hc = colorbar('eastoutside');
    m_grid('box','fancy'); m_coast('patch',nicecolor('wwk'));
    title('Percent change due to MLD_{max} (\deltaPOC flux/\deltaz x \deltaMLD_{max})')
subplot(325)
    m_proj('miller','lat',[latminplot latmaxplot],'lon',[lonminplot lonmaxplot])
    m_contourf([glon glon + 360], glat, repmat(100*(TaylorGrid.Residual_MLDmax_wMLDchange./NormalizeData),1,2), [cmin: cint: cmax],'linecolor','none'); hold on;
    colormap(C); caxis([cmin cmax]); hc = colorbar('eastoutside');
    m_grid('box','fancy'); m_coast('patch',nicecolor('wwk'));
    title('Residual of Taylor decomposition for change in MLD_{max} POC flux')
    
    %% Taylor decomposition at 100 m
figure(10); clf;
set(gcf,'color','w')
x0=0;
y0=0;
width=34;
height=14;
set(gcf,'units','centimeters','position',[x0,y0,width,height])
    %cmin = -75; cmax = 75;
    cmin = -2; cmax = 2; cint = 0.02;
    C = cmocean('balance');
    
    NormalizeData = 100; %outputGrid(1).POCflux_depthCriteria_AnnMean_clim_grid(:,:,5); %Annual mean POC flux from beginning of century

subplot(221)
    m_proj('miller','lat',[latminplot latmaxplot],'lon',[lonminplot lonmaxplot])
    m_contourf([glon glon + 360], glat, repmat(100*(TaylorGrid.dEPdt_100m./NormalizeData),1,2), [cmin: cint: cmax],'linecolor','none'); hold on;
    colormap(C); caxis([cmin cmax]); hc = colorbar('eastoutside'); ylabel(hc,'mol C m^{-2} yr^{-1}')
    m_grid('box','fancy'); m_coast('patch',nicecolor('wwk'));
    title('a) Change in POC flux evaluated at 100 m')
subplot(222)
    m_proj('miller','lat',[latminplot latmaxplot],'lon',[lonminplot lonmaxplot])
    m_contourf([glon glon + 360], glat, repmat(100*(TaylorGrid.NPPterm./NormalizeData),1,2), [cmin: cint: cmax],'linecolor','none'); hold on;
    colormap(C); caxis([cmin cmax]); hc = colorbar('eastoutside'); ylabel(hc,'mol C m^{-2} yr^{-1}')
    m_grid('box','fancy'); m_coast('patch',nicecolor('wwk'));
    title('b) NPP term evaluated at 100 m (\deltaNPP/\deltat x e-ratio)')
subplot(223)
    m_proj('miller','lat',[latminplot latmaxplot],'lon',[lonminplot lonmaxplot])
    m_contourf([glon glon + 360], glat, repmat(100*(TaylorGrid.eRatioTerm_100m./NormalizeData),1,2), [cmin: cint: cmax],'linecolor','none'); hold on;
    colormap(C); caxis([cmin cmax]); hc = colorbar('eastoutside'); ylabel(hc,'mol C m^{-2} yr^{-1}')
    m_grid('box','fancy'); m_coast('patch',nicecolor('wwk'));
    title('c) e-ratio term evaluated at (\deltae-ratio/\deltat x NPP)')
subplot(224)
    m_proj('miller','lat',[latminplot latmaxplot],'lon',[lonminplot lonmaxplot])
    m_contourf([glon glon + 360], glat, repmat(100*(TaylorGrid.Residual_100m./NormalizeData),1,2), [cmin: cint: cmax],'linecolor','none'); hold on;
    colormap(C); caxis([cmin cmax]); hc = colorbar('eastoutside'); ylabel(hc,'mol C m^{-2} yr^{-1}')
    m_grid('box','fancy'); m_coast('patch',nicecolor('wwk'));
    title('d) Residual evaluated at 100 m')
    
    %%
figure(11); clf
set(gcf,'color','w')
x0=0;
y0=0;
width=15;
height=14;
set(gcf,'units','centimeters','position',[x0,y0,width,height])
cmin = -1; cmax = 1; cint = 0.1;

subplot(2,1,1)
    %Ratio is normalized by absolute value of percent change
    ratio_plot = (TaylorGrid.eRatioTerm_100m./TaylorGrid.NPPterm).*abs(TaylorGrid.dEPdt_100m./outputGrid(1).POCflux_depthCriteria_AnnMean_clim_grid(:,:,5));
    [len, ~] = size(ratio_plot);
    for i = 1:len
        ind = find(ratio_plot(i,:) < cmin);
        ratio_plot(i,ind) = cmin;
        ind = find(ratio_plot(i,:) > cmax);
        ratio_plot(i,ind) = cmax;
    end
    m_proj('miller','lat',[latminplot latmaxplot],'lon',[lonminplot lonmaxplot])
    m_contourf([glon glon + 360], glat, repmat(ratio_plot,1,2), [cmin: cint: cmax],'linecolor','none'); hold on;
    colormap(C); caxis([cmin cmax]); hc = colorbar('eastoutside');
    m_grid('box','fancy'); m_coast('patch',nicecolor('wwk'));
    title('Ratio of NPP:e-ratio influence, 100 m')
subplot(2,1,2)
    ratio_plot = (TaylorGrid.eRatioTerm_MLDmax./TaylorGrid.NPPterm).*abs(TaylorGrid.dEPdt_MLDmax./outputGrid(1).POCflux_depthCriteria_AnnMean_clim_grid(:,:,4));
    [len, ~] = size(ratio_plot);
    for i = 1:len
        ind = find(ratio_plot(i,:) < cmin);
        ratio_plot(i,ind) = cmin;
        ind = find(ratio_plot(i,:) > cmax);
        ratio_plot(i,ind) = cmax;
    end
    m_proj('miller','lat',[latminplot latmaxplot],'lon',[lonminplot lonmaxplot])
    m_contourf([glon glon + 360], glat, repmat(ratio_plot,1,2), [cmin: cint: cmax],'linecolor','none'); hold on;
    colormap(C); caxis([cmin cmax]); hc = colorbar('eastoutside');
    m_grid('box','fancy'); m_coast('patch',nicecolor('wwk'));
    title('Ratio of NPP:e-ratio influence, MLD_{max}')
    
%% Plot baseline values for POC flux and e-ratio    
figure(12); clf
set(gcf,'color','w')
x0=0;
y0=0;
width=32;
height=14;
set(gcf,'units','centimeters','position',[x0,y0,width,height])
C = cmocean('deep');

subplot(2,2,1)
    cmin = 0; cmax = 6; cint = 0.2;
    m_proj('miller','lat',[latminplot latmaxplot],'lon',[lonminplot lonmaxplot])
    m_contourf([glon glon + 360], glat, repmat(outputGrid(1).POCflux_depthCriteria_AnnMean_clim_grid(:,:,5),1,2), [cmin: cint: cmax],'linecolor','none'); hold on;
    colormap(C); caxis([cmin cmax]); hc = colorbar('eastoutside');
    m_grid('box','fancy'); m_coast('patch',nicecolor('wwk'));
    title('Annual mean POC flux at 100 m, 2005-2024')
subplot(2,2,3)
    cmin = 0; cmax = 0.4; cint = 0.02;
    eratio_plot = (outputGrid(1).POCflux_depthCriteria_AnnMean_clim_grid(:,:,5)./outputGrid(1).NPPsum_depthCriteria_AnnMean_clim_grid(:,:,5));
    m_proj('miller','lat',[latminplot latmaxplot],'lon',[lonminplot lonmaxplot])
    m_contourf([glon glon + 360], glat, repmat(eratio_plot,1,2), [cmin: cint: cmax], 'linecolor','none'); hold on;
    colormap(C); caxis([cmin cmax]); hc = colorbar('eastoutside');
    m_grid('box','fancy'); m_coast('patch',nicecolor('wwk'));
    title('Annual mean e-ratio at 100 m, 2005-2024')
subplot(2,2,2)
    cmin = 0; cmax = 6; cint = 0.2;
    m_proj('miller','lat',[latminplot latmaxplot],'lon',[lonminplot lonmaxplot])
    m_contourf([glon glon + 360], glat, repmat(outputGrid(1).POCflux_depthCriteria_AnnMean_clim_grid(:,:,4),1,2), [cmin: cint: cmax],'linecolor','none'); hold on;
    colormap(C); caxis([cmin cmax]); hc = colorbar('eastoutside');
    m_grid('box','fancy'); m_coast('patch',nicecolor('wwk'));
    title('Annual mean POC flux at MLD_{max}, 2005-2024')
subplot(2,2,4)
    cmin = 0; cmax = 0.4; cint = 0.02;
    eratio_plot = (outputGrid(1).POCflux_depthCriteria_AnnMean_clim_grid(:,:,4)./outputGrid(1).NPPsum_depthCriteria_AnnMean_clim_grid(:,:,4));
    m_proj('miller','lat',[latminplot latmaxplot],'lon',[lonminplot lonmaxplot])
    m_contourf([glon glon + 360], glat, repmat(eratio_plot,1,2), [cmin: cint: cmax], 'linecolor','none'); hold on;
    colormap(C); caxis([cmin cmax]); hc = colorbar('eastoutside');
    m_grid('box','fancy'); m_coast('patch',nicecolor('wwk'));
    title('Annual mean e-ratio at MLD_{max}, 2005-2024')
    %%
figure(13); clf
set(gcf,'color','w')
x0=0;
y0=0;
width=16;
height=14;
set(gcf,'units','centimeters','position',[x0,y0,width,height])
C = cmocean('balance');

subplot(2,1,1)
    cmin = -1.2; cmax = 1.2; cint = 0.1;
    m_proj('miller','lat',[latminplot latmaxplot],'lon',[lonminplot lonmaxplot])
    m_contourf([glon glon + 360], glat, repmat(outputGrid(1).POCflux_depthCriteria_AnnMean_clim_grid(:,:,5) - outputGrid(1).POCflux_depthCriteria_AnnMean_clim_grid(:,:,4),1,2), [cmin: cint: cmax],'linecolor','none'); hold on;
    colormap(C); caxis([cmin cmax]); hc = colorbar('eastoutside');
    m_grid('box','fancy'); m_coast('patch',nicecolor('wwk'));
    title('POC flux difference, 100 m - MLD_{max}, 2005-2024')
subplot(2,1,2)
    cmin = -0.2; cmax = 0.2; cint = 0.02;
    eratio_plot = (outputGrid(1).POCflux_depthCriteria_AnnMean_clim_grid(:,:,5)./outputGrid(1).NPPsum_depthCriteria_AnnMean_clim_grid(:,:,5)) - ...
        (outputGrid(1).POCflux_depthCriteria_AnnMean_clim_grid(:,:,4)./outputGrid(1).NPPsum_depthCriteria_AnnMean_clim_grid(:,:,4));
    m_proj('miller','lat',[latminplot latmaxplot],'lon',[lonminplot lonmaxplot])
    m_contourf([glon glon + 360], glat, repmat(eratio_plot,1,2), [cmin: cint: cmax], 'linecolor','none'); hold on;
    colormap(C); caxis([cmin cmax]); hc = colorbar('eastoutside');
    m_grid('box','fancy'); m_coast('patch',nicecolor('wwk'));
    title('e-ratio difference, 100 m - MLD_{max}, 2005-2024')