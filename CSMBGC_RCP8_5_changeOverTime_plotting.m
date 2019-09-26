%% Regrid for plotting
%Note that if the regridded grid size gets too fine for the input data,
%this function creates anomalous patterns where input data is scarce
    %Also note that it might be better to make all calculations (e-ratio,
    %etc. prior to regridding)
lonmin = 0; lonmax = 360;
latmin = -85; latmax = 85;
lon_gridsize = 1; lat_gridsize = 1;

% Pack depth criteria, POC flux at depth criteria, and NPP sum at depth criteria into stack
for i = 1:numCriteria
    for j = 1:2
        stackedCriteria(:,:,i + 20*(j-1)) = output(j).depthCriteria_wAnnMean_clim{i}; %mld_all, zcomp, zeu, mldmax, 100 m
        stackedCriteria(:,:,i+numCriteria + 20*(j-1)) = output(j).POCflux_depthCriteria_AnnMean_clim{i};
        stackedCriteria(:,:,i+2*numCriteria + 20*(j-1)) = output(j).NPPsum_depthCriteria_AnnMean_clim{i};
        stackedCriteria(:,:,i+3*numCriteria + 20*(j-1)) = output(j).NPPdiat_depthCriteria_AnnMean_clim{i};
    end
end

% Add Taylor decomposition into stacked criteria
stackedCriteria(:,:,41) = Taylor_NPPTerm;
stackedCriteria(:,:,42) = Taylor_eRatioTerm;
stackedCriteria(:,:,43) = Taylor_eRatioTerm_MLDmax;
stackedCriteria(:,:,44) = Taylor_MLDchange;
stackedCriteria(:,:,45) = Taylor_Residual_100m;
stackedCriteria(:,:,46) = Taylor_Residual_MLDmax;
stackedCriteria(:,:,47) = dEPdt;
stackedCriteria(:,:,48) = dEPdt_MLDmax;
stackedCriteria(:,:,49) = Taylor_Residual_MLDmax_wMLDchange;

[glon, glat, stackedCriteria_grid] = regrid_even(TLONG, TLAT, stackedCriteria,...
    lonmin, lonmax, latmin, latmax, lon_gridsize, lat_gridsize);

%% Put back in output format
outputGridToday.depthCriteria_wAnnMean_clim_grid = stackedCriteria_grid(:,:,1:5);
outputGridToday.POCflux_depthCriteria_AnnMean_clim_grid = stackedCriteria_grid(:,:,6:10);
outputGridToday.NPPsum_depthCriteria_AnnMean_clim_grid = stackedCriteria_grid(:,:,11:15);
outputGridToday.NPPdiat_depthCriteria_AnnMean_clim_grid = stackedCriteria_grid(:,:,16:20);

outputGridEndCentury.depthCriteria_wAnnMean_clim_grid = stackedCriteria_grid(:,:,21:25);
outputGridEndCentury.POCflux_depthCriteria_AnnMean_clim_grid = stackedCriteria_grid(:,:,26:30);
outputGridEndCentury.NPPsum_depthCriteria_AnnMean_clim_grid = stackedCriteria_grid(:,:,31:35);
outputGridEndCentury.NPPdiat_depthCriteria_AnnMean_clim_grid = stackedCriteria_grid(:,:,36:40);

outputGrid = [outputGridToday, outputGridEndCentury];

TaylorGrid.NPPterm = stackedCriteria_grid(:,:,41);
TaylorGrid.eRatioTerm_100m = stackedCriteria_grid(:,:,42);
TaylorGrid.eRatioTerm_MLDmax = stackedCriteria_grid(:,:,43);
TaylorGrid.MLDchange = stackedCriteria_grid(:,:,44);
TaylorGrid.Residual_100m = stackedCriteria_grid(:,:,45);
TaylorGrid.Residual_MLDmax = stackedCriteria_grid(:,:,46);
TaylorGrid.dEPdt_100m = stackedCriteria_grid(:,:,47);
TaylorGrid.dEPdt_MLDmax = stackedCriteria_grid(:,:,48);
TaylorGrid.Residual_MLDmax_wMLDchange = stackedCriteria_grid(:,:,49);

%% Compare MLD_max for beginning and end of century

F = 10; %fontsize

%Maximum annual MLD at beginning and end of century
    figure(1); clf;
set(gcf,'color','w')
x0=5;
y0=5;
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
        title('Climatological maximum annual MLD, 2005-2024 (m)')
    elseif i == 2
        title('Climatological maximum annual MLD, 2081-2100 (m)')
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
    m_plot(lonstn, latstn, 'm.','markersize',20); hold on;
    colormap(C); caxis([cmin cmax]); hc = colorbar('eastoutside');
    m_grid('box','fancy'); m_coast('patch',nicecolor('wwk'));
    title('Change in climatological maximum annual MLD over 21st century (m)')
    
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
    %m_plot(lonstn, latstn, 'm.','markersize',20); hold on;
    colormap(C); caxis([cmin cmax]); hc = colorbar('eastoutside');
    m_grid('box','fancy'); m_coast('patch',nicecolor('wwk'));
    title('Percent change in climatological maximum annual MLD over 21st century (%)')
    
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
        title('POC flux change over 21st century, evaluated at 100 m')
    elseif i == 2
        title('POC flux change over 21st century, evaluated at maximum annual MLD')
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
    title('Difference in POC flux change, maximum annual MLD - 100 m')
    
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
        title('e-ratio change over 21st century, evaluated at 100 m')
    elseif i == 2
        title('e-ratio change over 21st century, evaluated at maximum annual MLD')
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
    title('Difference in e-ratio change, maximum annual MLD - 100 m')
    
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

fractDiatToday = outputGridToday.NPPdiat_depthCriteria_AnnMean_clim_grid(:,:,j)./outputGridToday.NPPsum_depthCriteria_AnnMean_clim_grid(:,:,j);
fractDiatEndCent = outputGridEndCentury.NPPdiat_depthCriteria_AnnMean_clim_grid(:,:,j)./outputGridEndCentury.NPPsum_depthCriteria_AnnMean_clim_grid(:,:,j);
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

  %%  
figure(9); clf;
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
    
%%
figure(10); clf;
set(gcf,'color','w')
x0=0;
y0=0;
width=34;
height=14;
set(gcf,'units','centimeters','position',[x0,y0,width,height])
    cmin = -75; cmax = 75;
    C = cmocean('balance');
    
    NormalizeData = outputGrid(1).POCflux_depthCriteria_AnnMean_clim_grid(:,:,5); %Annual mean POC flux from beginning of century

subplot(221)
    m_proj('miller','lat',[latminplot latmaxplot],'lon',[lonminplot lonmaxplot])
    m_contourf([glon glon + 360], glat, repmat(100*(TaylorGrid.dEPdt_100m./NormalizeData),1,2), [cmin: cint: cmax],'linecolor','none'); hold on;
    colormap(C); caxis([cmin cmax]); hc = colorbar('eastoutside');
    m_grid('box','fancy'); m_coast('patch',nicecolor('wwk'));
    title('Percent change in POC flux at 100 m, from 2005-2024 to 2081-2100')
subplot(222)
    m_proj('miller','lat',[latminplot latmaxplot],'lon',[lonminplot lonmaxplot])
    m_contourf([glon glon + 360], glat, repmat(100*(TaylorGrid.NPPterm./NormalizeData),1,2), [cmin: cint: cmax],'linecolor','none'); hold on;
    colormap(C); caxis([cmin cmax]); hc = colorbar('eastoutside');
    m_grid('box','fancy'); m_coast('patch',nicecolor('wwk'));
    title('Percent change due to NPP (\deltaNPP/\deltat x e-ratio)')
subplot(223)
    m_proj('miller','lat',[latminplot latmaxplot],'lon',[lonminplot lonmaxplot])
    m_contourf([glon glon + 360], glat, repmat(100*(TaylorGrid.eRatioTerm_100m./NormalizeData),1,2), [cmin: cint: cmax],'linecolor','none'); hold on;
    colormap(C); caxis([cmin cmax]); hc = colorbar('eastoutside');
    m_grid('box','fancy'); m_coast('patch',nicecolor('wwk'));
    title('Percent change due to e-ratio (\deltae-ratio/\deltat x NPP)')
subplot(224)
    m_proj('miller','lat',[latminplot latmaxplot],'lon',[lonminplot lonmaxplot])
    m_contourf([glon glon + 360], glat, repmat(100*(TaylorGrid.Residual_100m./NormalizeData),1,2), [cmin: cint: cmax],'linecolor','none'); hold on;
    colormap(C); caxis([cmin cmax]); hc = colorbar('eastoutside');
    m_grid('box','fancy'); m_coast('patch',nicecolor('wwk'));
    title('Residual of Taylor decomposition for change in 100 m POC flux')
    
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