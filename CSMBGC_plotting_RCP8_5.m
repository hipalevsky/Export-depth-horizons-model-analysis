% This script is to be run in CSMBGC_exportanalysis, using interpolated and
% regridded output.

fignum = 1; %make counter for figure numbers

%Set location of sections at beginning so can be shown on maps
lonAtl = 330; lonPac = 180;
sectionTitle = {'Atlantic Ocean (30^oW)', 'Pacific Ocean (180^oW)'};
latminplot_sect = -70; latmaxplot_sect = 60.5;
    L = 1.5; %linewidth for map plots
    MS = 6; % marker size for stations on map plots
    F = 10; %fontsize

%% Plot difference between POC flux at different depth criteria in map view

parameterTitle = {'NPP', 'POC flux', 'e-ratio (POC flux/NPP)'};
criteriaTitle = {'100 meters', 'Maximum annual MLD'};
for j = 2 %1:3
    if j == 1
criteriaToPlot1 = stackedCriteria_grid(:,:,16); % NPP at 100m, in mol C m-2 yr-1 (zeu = 13, 100m = 16)
criteriaToPlot2 = stackedCriteria_grid(:,:,15); % NPP at maximum annual MLD, in mol C m-2 yr-1
cmin = -8; cmax = 8; cint = 0.5;
    elseif j == 2
criteriaToPlot1 = stackedCriteria_grid(:,:,7); % POC flux at 100m, in mol C m-2 yr-1 (mld = 7, zeu = 8, 100m = 11)
criteriaToPlot2 = stackedCriteria_grid(:,:,10); % POC flux at maximum annual MLD, in mol C m-2 yr-1
cmin = -1.2; cmax = 1.2; cint = 0.05;
    elseif j == 3
criteriaToPlot1 = stackedCriteria_grid(:,:,11)./stackedCriteria_grid(:,:,16); % e-ratio at 100m
criteriaToPlot2 = stackedCriteria_grid(:,:,10)./stackedCriteria_grid(:,:,15); % e-ratio at maximum annual MLD
cmin = -0.15; cmax = 0.15; cint = 0.01;
    end

    figure(fignum); clf; fignum = fignum + 1;  
set(gcf,'color','w')
x0=5;
y0=5;
width=16;
height=7;
set(gcf,'units','centimeters','position',[x0,y0,width,height])

latminplot = -75; latmaxplot = 75;
lonminplot = 0; lonmaxplot = 360;
%C = (colormap(lbmap(61,'BrownBlue')));
%divergemap2; C = diverge_map2;
C = cmocean('balance');
set(0,'defaultAxesFontSize',F)
m_proj('miller','lat',[latminplot latmaxplot],'lon',[lonminplot lonmaxplot])
m_contourf(glon, glat, criteriaToPlot1 - criteriaToPlot2,[cmin: cint: cmax],'linecolor','none'); hold on;
colormap(C); caxis([cmin cmax]); hc = colorbar('eastoutside');
m_grid('box','fancy'); m_coast('patch',nicecolor('wwk'));
if j <= 2
    xlabel(hc,'mol C m^{-2} yr^{-1}','Fontsize',F);
end
%Add lines showing sections
m_plot([lonAtl lonAtl],[latminplot_sect latmaxplot_sect],'k--','linewidth',L); hold on;
m_plot([lonPac lonPac],[latminplot_sect latmaxplot_sect],'k--','linewidth',L); hold on;
m_plot(stn_plot(:,2), stn_plot(:,1),'ko','markerfacecolor','y','markersize',MS); hold on;
title([parameterTitle{j} ' difference, ', criteriaTitle{1}, ' - ', criteriaTitle{2}]);

text(-2.3,1.3,'b','Fontsize',F+2,'fontweight','bold')

end

%% Plot depth criteria   
criteriaTitle = {'Euphotic depth', 'Maximum annual mixed layer depth (MLD)','Particle compensation depth'};
    criteriaNum = [3 5 4];
cmin = 0; cmax = 340; cint = 20;
latminplot = -75; latmaxplot = 75;
lonminplot = 0; lonmaxplot = 360;
% Plot depths of the two individual criteria
for i = 2
    figure(fignum); clf; fignum = fignum + 1;  
set(gcf,'color','w')
x0=5;
y0=5;
width=16;
height=7;
set(gcf,'units','centimeters','position',[x0,y0,width,height])

m_proj('miller','lat',[latminplot latmaxplot],'lon',[lonminplot lonmaxplot])
m_contourf(glon, glat, stackedCriteria_grid(:,:,criteriaNum(i)),[cmin: cint: cmax],'linecolor','none'); hold on;
colormap(cmocean('deep')); caxis([cmin cmax]); hc = colorbar('eastoutside');
m_grid('box','fancy'); m_coast('patch',nicecolor('wwk'));
%Add lines showing sections
m_plot([lonAtl lonAtl],[latminplot_sect latmaxplot_sect],'k--','linewidth',L); hold on;
m_plot([lonPac lonPac],[latminplot_sect latmaxplot_sect],'k--','linewidth',L); hold on;
m_plot(stn_plot(:,2), stn_plot(:,1),'ko','markerfacecolor','y','markersize',MS); hold on;
title(criteriaTitle{i}); xlabel(hc,'Depth (m)','Fontsize',F);

text(-2.3,1.3,'a','Fontsize',F+2,'fontweight','bold')
end

%% Plot depth difference between the two criteria
% cmin = -200; cmax = 200; cint = 20;
%     figure(fignum); clf; fignum = fignum + 1;       
% m_proj('miller','lat',[latminplot latmaxplot],'lon',[lonminplot lonmaxplot])
% m_contourf(glon, glat, stackedCriteria_grid(:,:,criteriaNum(2)) - stackedCriteria_grid(:,:,criteriaNum(1)),[cmin: cint: cmax],'linecolor','none'); hold on;
% colormap(C); caxis([cmin cmax]); hc = colorbar('southoutside');
% m_grid('box','fancy'); m_coast('patch',nicecolor('wwk'));
% title([criteriaTitle{2}, ' - ', criteriaTitle{1}]); xlabel(hc,'Depth (m)','Fontsize',F);

%% Plot POC flux and e-ratio at depth criteria
% criteriaTitle = {'Seasonal MLD', 'Euphotic depth', 'Particle compensation depth', 'Maximum annual MLD', '100 meters'};
% cmin = 0; cmax = 3; cint = 0.1;
% % Plot depths of the two individual criteria
% for i = 1:5
%     figure(fignum); clf; fignum = fignum + 1;       
% m_proj('miller','lat',[latminplot latmaxplot],'lon',[lonminplot lonmaxplot])
% m_contourf(glon, glat, stackedCriteria_grid(:,:,i + numCriteria + 1),[cmin: cint: cmax],'linecolor','none'); hold on;
% colormap(C); caxis([cmin cmax]); hc = colorbar('southoutside');
% m_grid('box','fancy'); m_coast('patch',nicecolor('wwk'));
% title(['POC flux at ', criteriaTitle{i}]); xlabel(hc,'mol C m^{-2} yr^{-1}','Fontsize',F);
% end

%% Plot correlation between POC flux at different depth criteria and the depth criteria themselves
% Note that this should be done with native grid, not regridded data for plotting

    figure(fignum); clf; fignum = fignum + 1;
set(gcf,'color','w')
x0=5;
y0=5;
width=12;
height=7;
set(gcf,'units','centimeters','position',[x0,y0,width,height])

criteriaToCompare1 = stackedCriteria(:,:,11); % POC flux at 100 m (11), in mol C m-2 yr-1
criteriaToCompare2 = stackedCriteria(:,:,10); % POC flux at maximum annual MLD, in mol C m-2 yr-1
depthCriteriaForCorrelation = stackedCriteria(:,:,5); % maximum annual MLD
depthCriteriaTitle = {'Maximum annual MLD (m)'};
    depthtol = 60; % don't include samples where depthCriteriaForCorrelation is less than this value

xmin = depthtol - 1; xmax = 500;
ymin = -2; ymax = 2.9;
colormap(parula);

A = criteriaToCompare1 - criteriaToCompare2; A = A(:);
B = depthCriteriaForCorrelation(:);
ind = find(B > depthtol & isnan(A + B) == 0);
[rho,df,rho_sig95] = correlate(B(ind),A(ind));

dscatter(B(ind),A(ind)); hold on;
plot([1:1000],0*ones(1000,1),'k--'); hold on;
axis([xmin xmax ymin ymax])
xlabel(depthCriteriaTitle,'Fontsize',F)
%ylabel('mol C m^{-2} yr^{-1}','Fontsize',F)
ylabel({'POC flux difference'; '100 m - Maximum annual MLD'; 'mol C m^{-2} yr^{-1}'},'Fontsize',F)

text(72,2.5,'c','Fontsize',F+2,'fontweight','bold')

%% Plot POC flux in sections
% Note that handling of minimum annual MLD (more relevant than weighted
% mean seasonal MLD) is clunky

%Prepare sections
    [~, glonid_Atl] = min(abs(glon - lonAtl));
    [~, glonid_Pac] = min(abs(glon - lonPac));
    sections = [glonid_Atl glonid_Pac];
    latplotid = find(glat >= latminplot_sect & glat <= latmaxplot_sect);
zminplot = 0; zmaxplot = 245;
    zplotid = find(z >= zminplot & z <= zmaxplot);

%Version with each parameter in different contour plot
% parameterTitle = {'NPP', 'POC flux', 'e-ratio (POC flux/NPP)'};
% for j = 2:3
%     if j == 1
% dataForSection = NPPsum_interp_AnnMean_clim_grid;
% cmin = 0; cmax = 30; cint = 1;
%     elseif j == 2
% dataForSection = POCflux_interp_AnnMean_clim_grid;
% cmin = 0; cmax = 4; cint = 0.25;
%     elseif j == 3
% dataForSection = POCflux_interp_AnnMean_clim_grid./NPPsum_interp_AnnMean_clim_grid;
% cmin = 0; cmax = 0.35; cint = 0.025;
%     end
% 
% C = (colormap(lbmap(61,'Blue')));
% lineC = [nicecolor('mw'); nicecolor('mwww'); nicecolor('ryy'); nicecolor('gcck'); nicecolor('k'); nicecolor('kww')];
% Lwid = 2;
% 
%     figure(fignum); clf; fignum = fignum + 1;
% for i = 1:length(sections)
%     subplot(1,2,i)
%         contourf(glat(latplotid), z(zplotid), flipud(rot90(squeeze(dataForSection(latplotid,sections(i),zplotid)))),...
%             [cmin: cint: cmax],'linecolor','none'); hold on;
%         colormap(C); hc = colorbar; set(gca, 'YDir', 'reverse');
%         for k = [1,3:6]
%             plot(glat(latplotid),squeeze(stackedCriteria_grid(latplotid,sections(i),k)),'-','color',lineC(k,:),'linewidth',Lwid); hold on;
%         end
%         title([sectionTitle{i} ' section of ' parameterTitle{j}]);
%         if j <= 2
%             ylabel(hc,'mol C m^{-2} d^{-1}');
%         end
%         ylabel('Depth (m)');
%         if i == 2
%             xlabel('Latitude');
%         end
%         caxis([cmin cmax])
%         %Label depth criteria
% end
% 
% end

%Version with POC flux in color and e-ratio in dashes on same plot
dataForSection = POCflux_interp_AnnMean_clim_grid;
cmin = 0; cmax = 4; cint = 0.25;

dataForSectionb = POCflux_interp_AnnMean_clim_grid./NPPsum_interp_AnnMean_clim_grid;
cminb = 0; cmaxb = 0.35; cintb = 0.025; Clab = [0.05 0.1 0.15 0.2]; 

    figure(fignum); clf; fignum = fignum + 1;
set(gcf,'color','w')
x0=5;
y0=5;
width=21;
height=11.25;
set(gcf,'units','centimeters','position',[x0,y0,width,height])

C = (colormap(lbmap(61,'Blue')));
lineC = [nicecolor('mw'); nicecolor('mwww'); nicecolor('gcck'); nicecolor('ryy'); nicecolor('k'); nicecolor('kww')];
Lwid = 2.5;
for i = 1:length(sections)
    subplot(2,5,[2:5] + 5*(i-1))
        contourf(glat(latplotid), z(zplotid), flipud(rot90(squeeze(dataForSection(latplotid,sections(i),zplotid)))),...
            [cmin: cint: cmax],'linecolor','none'); hold on;
        colormap(C); hc = colorbar; set(gca, 'YDir', 'reverse'); ylabel(hc,'POC flux (mol C m^{-2} yr^{-1})')
        [Con,h] = contour(glat(latplotid), z(zplotid), flipud(rot90(squeeze(dataForSectionb(latplotid,sections(i),zplotid)))),...
            [cminb: cintb: cmaxb],'--','linecolor',nicecolor('kw'),'linewidth',0.75); hold on;
        clabel(Con,h,Clab,'Fontsize',8)
        for k = [1,6,3,4,5]
            hleg(k) = plot(glat(latplotid),squeeze(stackedCriteria_grid(latplotid,sections(i),k)),'-','color',lineC(k,:),'linewidth',Lwid); hold on;
        end
        title([sectionTitle{i} ' section']);
        ylabel('Depth (m)'); caxis([cmin cmax])
        if i == 2
            xlabel('Latitude'); 
        end
        %legend(hleg([1,3,4,5,6]),'Minimum annual MLD','Euphotic depth','Particle compensation depth','Maximum annual MLD','100 meters','location','westoutside')
end

% % %legend
subplot(2,5,1)
axis([0 10 0 10])
set(gca,'Visible','off')
Top = 2;
Left = -10;
LineLen = 3;
GapRight = 1;
GapVert = 1.5;
LegendLabels = {'','Minimum annual MLD','Particle compensation depth','Euphotic depth','Maximum annual MLD','100 meters'};
FSLeg = 10;

for i = [2:6]
    text(Left + LineLen + GapRight,Top-GapVert*(i-1),LegendLabels(i),'Fontsize',FSLeg)
    %leg = plot([Left, Top - GapVert*(i-1)],'o','color','k','linewidth',l,'markerfacecolor','g','markersize',5); set(leg,'clipping','off');
    if i == 2
        l = line([Left Left + LineLen],[Top-GapVert*(i-1) Top-GapVert*(i-1)],'color',lineC(1,:),'linewidth',Lwid); set(l,'clipping','off');
    else
        l = line([Left Left + LineLen],[Top-GapVert*(i-1) Top-GapVert*(i-1)],'color',lineC(i,:),'linewidth',Lwid); set(l,'clipping','off');
    end
end

%% Plot along sections
%     figure(fignum); clf; fignum = fignum + 1; 
% lineC = [nicecolor('mw'); nicecolor('ryy'); nicecolor('gcck'); nicecolor('k'); nicecolor('kww')];
% Lwid = 2;
% for j = 1:2
%     for i = 1:length(sections)
%         subplot(2,2,i + 2*(j-1))
%             for k = [1,5,2,3,4]
%                 if j == 1
%                     hleg(k) = plot(glat(latplotid)',squeeze(stackedCriteria_grid(latplotid,sections(i),k + 6)),'-','color',lineC(k,:),'linewidth',Lwid); hold on;
%                     cmin = 0; cmax = 6; yTitle = {'POC flux'};
%                 elseif j == 2
%                     hleg(k) = plot(glat(latplotid)',squeeze(stackedCriteria_grid(latplotid,sections(i),k + 6))./squeeze(stackedCriteria_grid(latplotid,sections(i),k + 11))...
%                         ,'-','color',lineC(k,:),'linewidth',Lwid); hold on;
%                     cmin = 0; cmax = 0.3; yTitle = {'e-ratio'};
%                 end
%             end
%             if j == 1
%                 title([sectionTitle{i} ' section']);
%             end
%             xlabel('Latitude'); xlim([min(glat(latplotid)) max(glat(latplotid))]); ylim([cmin cmax])
%             ylabel(yTitle); 
%             %legend(hleg([1,5,2,3,4]),'Seasonal MLD','100 meters','Euphotic depth','Compensation depth','Maximum annual MLD','location','westoutside')
%     end
% end

%% Plot locations of discrete station points
%     figure(fignum); clf; fignum = fignum + 1;  
% latminplot = -75; latmaxplot = 75;
% lonminplot = 0; lonmaxplot = 360;
% cdiv = [1: 0.5: 3];
% C = colormap(jet); C = C(1:floor(length(C)/(length(cdiv))):end,:);
% F = 12;
% 
% m_proj('miller','lat',[latminplot latmaxplot],'lon',[lonminplot lonmaxplot])
% m_contourf(glon, glat, depthregID_grid,cdiv,'linecolor','none'); hold on;
% colormap(C(1:length(cdiv),:)); caxis([min(cdiv) max(cdiv)]); hc = colorbar('southoutside');
% m_grid('box','fancy'); m_coast('patch',nicecolor('wwk'));
% m_plot(stns(1:9,2),stns(1:9,1),'mp','markerfacecolor','m','markersize',15); hold on;

%% Pull out data at individual locations to plot seasonal data through year
%Should do this as a function where will put in station locations and data
%to pull time series from and out will come some plots.

% Find nearest point in CSM grid for each of the time series stations
for i = 1:length(stns)
    [inds(i,:),gridpt(i,:)] = findNearestPoint(TLAT,TLONG,stns(i,1),stns(i,2));
end

% Plot CSM time series at each station location
tspan = [1:16];
tlabels = {'J','F','M','A','M','J','J','A','S','O','N','D','J','F','M','A',};
C = (colormap(lbmap(61,'Blue')));
lineC = [nicecolor('mw'); nicecolor('gcck'); nicecolor('ryy'); nicecolor('k'); nicecolor('kww')];
Lwid = 2;
cmin = zeros(1,length(stns));
cmax = [5 0.1 15 0.3];
  
stackedForTSPlotting{1} = POCprod_all;
stackedForTSPlotting{2} = POCflux_all;
parameterTitle = {'POC production, mol C m^{-3} yr^{-1}', 'POC remineralization, mol C m^{-3} yr^{-1}', 'POC flux, mol C m^{-2} yr^{-1}', 'e-ratio (POC flux/NPP)'}; %'NPP (depth-integrated), mol C m^{-2} yr^{-1}', 
    numpars = length(stackedForTSPlotting);

for i = 1:length(stns)
%     figure(20 + i); clf
%     set(gcf,'color','w')
    x0=1;
    y0=1;
    width=22;
    height=20;
    set(gcf,'units','centimeters','position',[x0,y0,width,height])
zminplot = 0; zmaxplot = zmaxstn(i); ztickint = 100;
    zplotid = find(z >= zminplot & z <= zmaxplot);
    zint = z(2) - z(1);
    ztick = [min(zplotid)+1:ztickint/zint:max(zplotid)];
    zticklabel = z(ztick);
for k = 1:numpars
%subplot(numpars,1,k)
%imagesc(squeeze(stackedForTSPlotting{k}(inds(i,1),inds(i,2),zplotid,tspan))); hold on;
%colormap(C); hc = colorbar; title([stnNames(i) parameterTitle(k)]); xticks(tspan); xticklabels(tlabels); yticks(ztick); yticklabels(zticklabel);
    for j = 1:5
        %plot(squeeze(depthCriteria{j}(inds(i,1),inds(i,2),tspan))/zint,'-','color',lineC(j,:),'linewidth',Lwid); hold on;
        if j > 1 & j < 4
            %plot(tspan,repmat(squeeze(depthCriteria_wAnnMean{j}(inds(i,1),inds(i,2),ceil(tspan(1)/12)))/zint,length(tspan),1),'--','color',lineC(j,:),'linewidth',Lwid); hold on;
        end
        if k == 1;
            NPP_vals(i,j) = NPPsum_depthCriteria_AnnMean_clim{j}(inds(i,1),inds(i,2));
            POCflux_vals(i,j) = POCflux_depthCriteria_AnnMean_clim{j}(inds(i,1),inds(i,2));
        end
    end
    ylabel('Depth (m)'); set(gca,'ydir','reverse'); xlim([min(tspan) max(tspan)]);
    %caxis([cmin(k) cmax(k)]);
end

end

%% Plot the spread among the different depth criteria at the stations of interest
    cmToMeters = 1/100;
    molToPg = 12*10^-15;
    c = 0; %counter
figure(fignum); clf; fignum = fignum + 1;  
for k = 2:3
    subplot(1,3,k-1)
if k == 1;
    plotting = NPP_vals; labeling = {'Annual NPP, mol C m^{-2} yr^{-1}'};
elseif k == 2;
    plotting = [POCflux_sum_glob/(molToPg)/(nansum(nansum(TAREA*(cmToMeters)^2))); POCflux_vals]; labeling = {'Annual POC flux'}; units = {'mol C m^{-2} yr^{-1}'}; ymax = 4;
elseif k == 3;
    plotting = [POCflux_sum_glob./NPP_sum_glob; POCflux_vals./NPP_vals]; labeling = {'Annual e-ratio'}; units = {'POC flux/NPP'}; ymax = 0.22;
end
    for i = [1,5,2,3,4]
        c = c+1;
        plot([1:5],plotting([1,10,6,3,7],i),'.','color',lineC(i,:),'markersize',40-c); hold on;
    end
    xlim([0.6 5.4]); xticks([1:5]); xticklabels({'Global','SO','Kuroshio','Irminger','EqPac'}); title(labeling); ylabel(units);
    ylim([0 ymax]);
    if k == 3
        legend('Seasonal MLD','100 meters','Compensation depth','Euphotic depth','Maximum annual MLD')
    end
end

%% Plot zonal mean values and global values
cmToMeters = 1/100;
    figure(fignum); clf; fignum = fignum + 1;  
set(gcf,'color','w')
x0=5;
y0=5;
width=20;
height=11.25;
set(gcf,'units','centimeters','position',[x0,y0,width,height])

    molToPg = 12*10^-15;
    buffer = 0.1; %for plotting
    int = 3; % for plotting
    lineC = [nicecolor('mw'); nicecolor('gcck'); nicecolor('ryy'); nicecolor('k'); nicecolor('kww')];
    Lwid = 4.5-2;
    FS = 10;
    stnids = [10,6,3,7];
    %Get ids to plot stns on zonal mean grid
    for i = 1:length(stnids)
        [~,ind] = min(abs(latdiv - stns(stnids(i)-1,1)));
        a = latdiv(ind) - stns(stnids(i)-1,1);
        stnid_plot(i) = ind - a/latint;
    end
for k = 2:3
if k == 2
    plotting = POCflux_sum_reg./(repmat(area_reg,5,1)')/molToPg;    %/repmat(POCflux_sum_glob,7,1); %mean annual POC flux, mol C m-2 yr-1
    plotting_stns = [POCflux_sum_glob/(molToPg)/(nansum(nansum(TAREA*(cmToMeters)^2))); POCflux_vals]; 
    titleText = {'Zonal mean POC flux'};
    labelText = {'mol C m^{-2} yr^{-1}'};
    ymin = 0; ymax = 4.8;
    plotnum = [1:3];
elseif k == 3
    plotting = POCflux_sum_reg./NPP_sum_reg; %mean annual e-ratio
    plotting_stns = [POCflux_sum_glob./NPP_sum_glob; POCflux_vals./NPP_vals];
    titleText = {'Zonal mean e-ratio'};
    labelText = {'POC flux/NPP'};
    ymin = 0; ymax = 0.4;
    plotnum = [5:7];
end
    subplot(2,4,plotnum)
    for i = [1,5,2,3,4]
        plot([1:length(area_reg)],plotting(:,i),'-','color',lineC(i,:),'linewidth',Lwid); hold on;
    end
    for i = [1,5,2,3,4]
        plot(stnid_plot,plotting_stns(stnids,i),'o','markerfacecolor',lineC(i,:),'markeredgecolor','k','markersize',8); hold on;
    end
    xlim([1 - buffer length(area_reg) + buffer]); xticks([2:int:length(area_reg)-1] + 0.5);
    xticklabels(latdiv([2:int:length(area_reg)-1]) + latint/2); ylim([ymin ymax]);
    ylabel(labelText,'Fontsize',FS);
    title(titleText,'Fontsize',FS);
    set(gca,'Fontsize',FS);
    if k == 3
        xlabel('Latitude','Fontsize',FS)
        %legend('Seasonal MLD','100 meters','Euphotic depth','Compensation depth','Maximum annual MLD')
    end
end

%legend
Top = 0.62;
Left = length(area_reg) + 1.2;
LineLen = 1.8;
GapRight = 0.4;
GapVert = 0.04;
LegendLabels = {'Seasonal MLD','Particle compensation depth','Euphotic depth','Maximum annual MLD','100 meters'};
FSLeg = 10;

for i = 1:5
    l = line([Left Left + LineLen],[Top-GapVert*(i-1) Top-GapVert*(i-1)],'color',lineC(i,:),'linewidth',Lwid); set(l,'clipping','off');
    text(Left + LineLen + GapRight,Top-GapVert*(i-1),LegendLabels(i),'Fontsize',FSLeg)
end
    


