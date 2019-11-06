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
stackedCriteria(:,:,50) = Taylor_NPPTerm_MLDmaxBegCent;
stackedCriteria(:,:,51) = Taylor_eRatioTerm_MLDmaxBegCent;
stackedCriteria(:,:,52) = Taylor_dEPdt_MLDmaxChange;
stackedCriteria(:,:,53) = Taylor_Residual_MLDmax_BegVsEndCentCriteria;

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
TaylorGrid.NPPTerm_MLDmaxBegCent = stackedCriteria_grid(:,:,50);
TaylorGrid.eRatioTerm_MLDmaxBegCent = stackedCriteria_grid(:,:,51);
TaylorGrid.dEPdt_MLDmaxChange = stackedCriteria_grid(:,:,52);
TaylorGrid.Residual_MLDmax_BegVsEndCentCriteria = stackedCriteria_grid(:,:,53);

%% Save gridded output for plotting
save RCP8_5_GriddedOutput outputGrid TaylorGrid glon glat