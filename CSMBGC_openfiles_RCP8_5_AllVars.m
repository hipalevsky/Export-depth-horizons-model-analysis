function [out] = CSMBGC_openfiles_RCP8_5_AllVars(year); 

%Modified from CSMBGC_openfiles_RCP8_5 to save output as structure and
%include more variables

filename = ['b40.rcp8_5.1deg.bdrd.001.' num2str(year) '.nc']; %ncdisp(filename)

%Properties that are the same for all years
        time = ncread(filename,'time'); %'days since 0000-01-01 00:00:00'
        z_t = ncread(filename,'z_t'); %'depth from surface to midpoint of layer', 'centimeters'
        z_w = ncread(filename,'z_w'); %'depth from surface to top of layer', 'centimeters'
        out.TLAT = ncread(filename,'TLAT'); %'degrees_north' (size is nlon x nlat)
        out.TLONG = ncread(filename,'TLONG'); %'degrees_east' (size is nlon x nlat)
        dz = ncread(filename,'dz'); %'thickness of layer k', 'centimeters'
        dzw = ncread(filename,'dzw'); %'midpoint of k to midpoint of k+1', 'centimeters'
        out.REGION_MASK = ncread(filename,'REGION_MASK'); %'basin index number (signed integers)' (size is nlon x nlat)
        out.TAREA = ncread(filename,'TAREA'); %'area of T cells','centimeter^2' (size is nlon x nlat)
        %time_bound = ncread(filename,'time_bound'); %'boundaries for time-averaging interval','days since 0000-01-01 00:00:00'    
%Convert units
        secinyr = 60*60*24*365; %set constant
        out.z = z_t/100; %depth from surface to midpoint of layer in meters
        out.z_top = z_w/100; %depth from surface to top of layer in meters
        
%Properties that vary - physical variables
    out.TEMP = ncread(filename,'TEMP'); %Potential Temperature in degC (size is nlon x nlat x z_t x time) 
    out.SALT = ncread(filename,'SALT'); %Salinity in gram/kilogram (size is nlon x nlat x z_t x time)  
    HMXL = ncread(filename,'HMXL'); %'Mixed Layer Depth in centimeters (size is nlon x nlat x time)   

%Properties that vary - NPP and biomass
    photoC_diat = ncread(filename,'photoC_diat'); %C uptake by diatoms in mmol C/m^3/sec (size is nlon x nlat x z_t(1:15) x time)   
    photoC_diaz = ncread(filename,'photoC_diaz'); %C uptake by diazotrophs in mmol C/m^3/sec (size is nlon x nlat x z_t(1:15) x time)   
    photoC_sp = ncread(filename,'photoC_sp'); %C uptake by small phyto in mmol C/m^3/sec (size is nlon x nlat x z_t(1:15) x time)      
%     diatC = ncread(filename,'diatC'); %Diatom Carbon in mmol C/m^3 (size is nlon x nlat x z_t(1:15) x time)   
%     diazC = ncread(filename,'diazC'); %Diazotroph Carbon in mmol C/m^3 (size is nlon x nlat x z_t(1:15) x time)        
%     spC = ncread(filename,'spC'); %Small Phytoplankton Carbon in mmol C/m^3 (size is nlon x nlat x z_t(1:15) x time)      

%Properties that vary - Chemical tracers and gas fluxes    
     out.ALK = ncread(filename,'ALK'); %Alkalinity in meq/m^3, (size is nlon x nlat x z_t x time)
     out.DIC = ncread(filename,'DIC'); %DIC in mmol C/m^3 (size is nlon x nlat x z_t x time)
%     DOC = ncread(filename,'DOC'); %Dissolved Organic Carbon in mmol C/m^3 (size is nlon x nlat x z_t x time) 
     out.pCO2SURF = ncread(filename,'pCO2SURF'); %Surface pCO2 in ppmv (size is nlon x nlat x time)  
     FG_CO2 = ncread(filename,'FG_CO2'); %Carbon Dioxide Flux in mmol DIC/m^3 cm/sec (size is nlon x nlat x time)
     FG_O2 = ncread(filename,'STF_O2'); %'Oxygen Flux in mmol O2/m^3 cm/sec (size is nlon x nlat x time)
     out.ATM_CO2 = ncread(filename,'ATM_CO2'); %Atmospheric CO2 in ppmv (size is nlon x nlat x time) %NEW - not in Lima 2014 output
     out.O2 = ncread(filename,'O2'); %Dissolved oxygen in mmol/m^3 (size is nlon x nlat x z_t x time) %NEW - not in Lima 2014 output
    
%Properties that vary - POC flux
    POC_FLUX_IN = ncread(filename,'POC_FLUX_IN'); %Incoming Flux of POC in mmol POC/m^3 cm/sec (size is nlon x nlat x z_t x time)   
    POC_PROD = ncread(filename,'POC_PROD'); %Production of POC in mmol POC/m^3/sec (size is nlon x nlat x z_t x time)    
    
%Convert units for properties that vary
%     year = floor((time-1)/365); %the "-1" keeps last day of Dec in current yr
%         yrslist = unique(year);
%         day = (time-1) - year*365 + 1;
    out.mld = HMXL/100; %mixed layer depth in meters

    % POC data (flux, production)
    out.POCflux = POC_FLUX_IN.*secinyr/100/1000; %mmol POC/m^3 cm/sec to mol C m-2 yr-1
    out.POCprod = POC_PROD.*secinyr/1000; %mmol POC/m^3/sec to mol C m-3 yr-1
        [d1s,d2s,d3s,d4s] = size(out.POCflux); %size of dimensions

    % NPP data (diatom, diazotroph, and small phytoplankton groups)
        % Volumetric rates
    NPPv_diat = photoC_diat.*secinyr/1000; %mmol C/m^3/sec to mol C m-3 yr-1
    NPPv_diaz = photoC_diaz.*secinyr/1000; %mmol C/m^3/sec to mol C m-3 yr-1
    NPPv_sp = photoC_sp.*secinyr/1000; %mmol C/m^3/sec to mol C m-3 yr-1
        % Depth-integrated rates (note: calculated before interpolation to
        % avoid artefacts with NaNs in surface interpolation)
        dz_grid = shiftdim(repmat(dz(1:15)/100,1,d4s,d1s,d2s),2); %create dz grid the size of the NPP grid and convert units from cm to m
    out.NPP_diat = cumsum(NPPv_diat.*dz_grid, 3); %mol C m-2 yr-1
    out.NPP_diaz = cumsum(NPPv_diaz.*dz_grid, 3); %mol C m-2 yr-1
    out.NPP_sp = cumsum(NPPv_sp.*dz_grid, 3); %mol C m-2 yr-1
    out.NPP_sum = out.NPP_diat + out.NPP_diaz + out.NPP_sp;
    
     % Air-sea fluxesclea
     out.F_O2 = FG_O2.*secinyr/100/1000; %mmol O2/m^3 cm/sec to mol O2 m-2 yr-1
     out.F_CO2 = FG_CO2.*secinyr/100/1000; %mmol CO2/m^3 cm/sec to mol CO2 m-2 yr-1
    
end

   

