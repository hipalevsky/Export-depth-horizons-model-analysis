%Compare Lima et al. 2014 output with observations in response to reviewer
%request for GRL paper

% First do the main paper analysis to pull out CSMBGC data
%CSMBGC_exportanalysis

%% Read in observations of POC flux for comparison
%Marsay et al. 2015, PNAS
[POCFluxObsData,~,~] = xlsread('POCfluxData_forGRLModelValidation.xlsx');

%% Pull out POC flux profiles from model to go with each obs station
stnid = POCFluxObsData(:,1);
figure(200); clf
    set(0,'defaultAxesFontSize',12)
    M = 18; L = 1.5;
    C_stns = [nicecolor('rrrrk'); nicecolor('ry'); nicecolor('ryyw');...
        nicecolor('gggyk'); [0 0 0]; nicecolor('gcckw'); nicecolor('ccb');...
        nicecolor('bbbk'); nicecolor('br'); nicecolor('mrw')];
    titles = {'PAP','Irminger','Iceland','NSAG','','K2','ALOHA','EQPAC','KIWI','OSP'};
for i = [1:4,6:10]
    ids = find(stnid == i);
        latstn(i) = POCFluxObsData(min(ids),3);
        lonstn(i) = POCFluxObsData(min(ids),4);
        month = POCFluxObsData(min(ids),2);
    % Find nearest point in CCSM grid for each of the time series stations
    [inds(i,:),gridpt(i,:)] = findNearestPoint(TLAT,TLONG,latstn(i),lonstn(i));
    % Pull out POC flux profile from CCSM
    POCflux_modelprofile(:,i) = squeeze(POCflux_interp_clim(inds(i,1), inds(i,2), :, month));
    % Plot POC flux from model and from obs
    subplot(2,5,i)
    plot(POCflux_modelprofile(:,i)*1000/365,z_new,'k-','linewidth',L); hold on;
    errorbar(POCFluxObsData(ids,6),POCFluxObsData(ids,5),POCFluxObsData(ids,7),'Horizontal','.','color',C_stns(i,:),'markersize',M); hold on;
    if i == 1 | i == 6
        ylabel('Depth (m)');
    end
    set(gca,'ydir','reverse'); ylim([0 800])
    title(titles(i))
    %xlabel('POC flux (mmol C m^{-2} d^{-1})')
end

%% Plot a map of the stations
figure(1); clf
set(0,'defaultAxesFontSize',8)
m_proj('miller','lat',[-70 70],'lon',[90 360])
for i = [1:4,6:10];
    m_plot(lonstn(i),latstn(i),'.','color',C_stns(i,:),'markersize',M); hold on;
end
m_grid('box','fancy'); m_coast('patch',nicecolor('wwk'));

%% Next steps: Prettifying figure
% Make a map that shows all plotted stations --> try using colored dots to
% represent each station and then use those colors in the profiles

%% Use nlinfit to calculate z* (remin length scale)
% These are all very close together, so this doesn't seem like a helpful
% way to compare model with obs

% for i = [1:5,9:max(stnid)]
%     clear Fz F0 z_fit
%     % Fit from particle compensation depth down
%         [profilemax, maxind]  = max(POCflux_modelprofile(:,i));
%     %Specify the values to use to fit the remin length scale
%             bottom = 120; %only use profile down to 600 m
%         Fz = POCflux_modelprofile(maxind:bottom,i);
%         F0 = profilemax;
%         z_fit = z_new(maxind:bottom);
%     %Remin length scale function (see eqn 2 Boyd and Buesseler 2009, among others --> BB09 is missing negative)
%         remin_fun = @(zstar, z_fit) (F0.*exp(-(z_fit - z_fit(1))./zstar));
%     %Specify initial guess for fitting the function
%         initial = 500;
%     %Calculate z* for the model profile
%         zstar_model(i) = nlinfit(z_fit', Fz, remin_fun, initial);
%         err = rms(Fz - remin_fun(zstar_model(i), z_fit));
%     %Plot resulting fit and original model profile
%         figure(200 + i); clf
%         plot(remin_fun(zstar_model(i), z_fit), z_fit,'k-'); hold on;
%         plot(Fz, z_fit,'b.'); hold on;
%         ylabel('Depth (m)'); set(gca,'ydir','reverse');
% end
