%% Visualize the weird CSM grid
%Need to already have TLONG and TLAT in workspace

%Make vectors with all CSM grid points
LONGLIST = TLONG(:);
LATLIST = TLAT(:);

%Make vectors with even grid points
glon = [0.5:2:359.5];
glat = [-89.5:2:89.5];
[longrid,latgrid] = meshgrid(glon,glat);
    latgrid = flipud(latgrid);
latgridlist = latgrid(:);
longridlist = longrid(:);

%% Plot globe on molleweide projection
figure(1); clf
set(gcf,'color','w')

x0=2;
y0=2;
width=28;
height=18;
set(gcf,'units','centimeters','position',[x0,y0,width,height])

m_proj('Mollweide','lat',[-85 85],'lon',[0 360])
m_plot(LONGLIST,LATLIST,'r.','markersize',5); hold on;
%m_plot(longridlist,latgridlist,'b.','markersize',5); hold on;
m_grid('box','fancy')
m_coast('patch','k');

%% Plot northern hemisphere in mercator projection
figure(2); clf
set(gcf,'color','w')

x0=2;
y0=2;
width=28;
height=18;
set(gcf,'units','centimeters','position',[x0,y0,width,height])

m_proj('mercator','lat',[1 76],'lon',[0 360])
m_plot(LONGLIST,LATLIST,'r.','markersize',5); hold on;
%m_plot(longridlist,latgridlist,'b.','markersize',5); hold on;
m_grid('box','fancy')
m_coast('patch','k');