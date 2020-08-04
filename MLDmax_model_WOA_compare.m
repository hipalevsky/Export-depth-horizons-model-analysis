%% Compare climatological MLDs from CSMBGC project model output w/ WOA data

%% Load World Ocean Atlas 2013 MLD climatology
% Calculated with WorldOceanAtlasProcess2, using 0.125 density cutoff to
% calculate MLD from T and S data
filename = 'C:/Users/palevsky/Dropbox/MATLAB/WorldOceanAtlasOutput.mat';
load(filename, '-mat', 'mld1deg')
mldmax_woa13 = max(mld1deg,[],3);

%% Compare with beginning century model output

F = 10; %fontsize

%Maximum annual MLD at beginning and end of century
    figure(20); clf;
set(gcf,'color','w')
x0=2;
y0=2;
width=16;
height=14;
set(gcf,'units','centimeters','position',[x0,y0,width,height])

latminplot = -75; latmaxplot = 75;
lonminplot = 10; lonmaxplot = 390;
cmin = 0; cmax = 500; cint = 20;
C = cmocean('deep');
set(0,'defaultAxesFontSize',F)

    subplot(2,1,1)
    m_proj('miller','lat',[latminplot latmaxplot],'lon',[lonminplot lonmaxplot])
    m_contourf([glon glon + 360], glat, repmat(outputGrid(1).depthCriteria_wAnnMean_clim_grid(:,:,4),1,2),[cmin: cint: cmax],'linecolor','none'); hold on;
    colormap(C); caxis([cmin cmax]); hc = colorbar('eastoutside');
    m_grid('box','fancy'); m_coast('patch',nicecolor('wwk'));
        title('CESM1-BEC climatological MLD_{max}, 2005-2024 (m)')
        
    subplot(2,1,2)
    m_proj('miller','lat',[latminplot latmaxplot],'lon',[lonminplot lonmaxplot])
    m_contourf([0.5:1:720], [-89.5:89.5], repmat(mldmax_woa13',1,2),[cmin: cint: cmax],'linecolor','none'); hold on;
    colormap(C); caxis([cmin cmax]); hc = colorbar('eastoutside');
    m_grid('box','fancy'); m_coast('patch',nicecolor('wwk'));
        title('WOA13 climatological MLD_{max} (m)')
%%        
    figure(21); clf;
set(gcf,'color','w')
x0=5;
y0=5;
width=16;
height=7;
set(gcf,'units','centimeters','position',[x0,y0,width,height])

cmin = -400; cmax = 400; cint = 10;
C = cmocean('balance');

    m_proj('miller','lat',[latminplot latmaxplot],'lon',[lonminplot lonmaxplot])
    m_contourf([glon(1:end-1) glon(1:end-1) + 360], glat, repmat([outputGrid(1).depthCriteria_wAnnMean_clim_grid(:,1:end-1,4) - mldmax_woa13(:,5:end-5)'],1,2),[cmin: cint: cmax],'linecolor','none'); hold on;
    colormap(C); caxis([cmin cmax]); hc = colorbar('eastoutside');
    m_grid('box','fancy'); m_coast('patch',nicecolor('wwk'));
        title('CESM1-BEC - WOA13 MLD_{max} (m)')
        
%% Compare beginning of century MLD_max whether > or < 100 m

F = 10; %fontsize

%Maximum annual MLD at beginning and end of century
    figure(22); clf;
set(gcf,'color','w')
x0=2;
y0=2;
width=16;
height=14;
set(gcf,'units','centimeters','position',[x0,y0,width,height])

latminplot = -75; latmaxplot = 75;
lonminplot = 10; lonmaxplot = 390;
cmin = -100; cmax = 100; cint = 5;
C = cmocean('balance');
set(0,'defaultAxesFontSize',F)

    m_proj('miller','lat',[latminplot latmaxplot],'lon',[lonminplot lonmaxplot])
    m_contourf([glon glon + 360], glat, repmat(outputGrid(1).depthCriteria_wAnnMean_clim_grid(:,:,4),1,2) - 100,[cmin: cint: cmax],'linecolor','none'); hold on;
    colormap(C); caxis([cmin cmax]); hc = colorbar('eastoutside');
    m_grid('box','fancy'); m_coast('patch',nicecolor('wwk'));
        title('CESM-BEC climatological MLD_{max}, 2005-2024, difference from 100 (m)')

