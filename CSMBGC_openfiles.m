%% Code to open monthly files from CSMBGC run (last 10 years of Lima et al 2014)
% Files extracted from ftp://ftp.whoi.edu/whoinet/CSMBGC/hilary/

%% All variables, grouped by type (so some can be omitted)

%%%% Basic physical variables %%%%
filename = ['SALT_monthly_0833-0842.nc']; ncdisp(filename)
    SALT = ncread(filename,'SALT'); %Salinity in gram/kilogram (size is nlon x nlat x z_t x time)  
        time = ncread(filename,'time'); %'days since 0000-01-01 00:00:00'
        z_t = ncread(filename,'z_t'); %'depth from surface to midpoint of layer', 'centimeters'
        z_w = ncread(filename,'z_w'); %'depth from surface to top of layer', 'centimeters'
        TLAT = ncread(filename,'TLAT'); %'degrees_north' (size is nlon x nlat)
        TLONG = ncread(filename,'TLONG'); %'degrees_east' (size is nlon x nlat)
        dz = ncread(filename,'dz'); %'thickness of layer k', 'centimeters'
        dzw = ncread(filename,'dzw'); %'midpoint of k to midpoint of k+1', 'centimeters'
        REGION_MASK = ncread(filename,'REGION_MASK'); %'basin index number (signed integers)' (size is nlon x nlat)
        TAREA = ncread(filename,'TAREA'); %'area of T cells','centimeter^2' (size is nlon x nlat)
        time_bound = ncread(filename,'time_bound'); %'boundaries for time-averaging interval','days since 0000-01-01 00:00:00'    
filename = ['TEMP_monthly_0833-0842.nc']; ncdisp(filename)
    TEMP = ncread(filename,'TEMP'); %Potential Temperature in degC (size is nlon x nlat x z_t x time)   
filename = ['HMXL_monthly_0833-0842.nc']; ncdisp(filename)
    HMXL = ncread(filename,'HMXL'); %'Mixed Layer Depth in centimeters (size is nlon x nlat x time)    
 
%%%% NPP and biomass (by phyto types) %%%%   
filename = ['photoC_diat_monthly_0833-0842.nc']; ncdisp(filename)
    photoC_diat = ncread(filename,'photoC_diat'); %C uptake by diatoms in mmol C/m^3/sec (size is nlon x nlat x z_t x time)   
filename = ['photoC_diaz_monthly_0833-0842.nc']; ncdisp(filename)
    photoC_diaz = ncread(filename,'photoC_diaz'); %C uptake by diazotrophs in mmol C/m^3/sec (size is nlon x nlat x z_t x time)   
filename = ['photoC_sp_monthly_0833-0842.nc']; ncdisp(filename)
    photoC_sp = ncread(filename,'photoC_sp'); %C uptake by small phyto in mmol C/m^3/sec (size is nlon x nlat x z_t x time)  
    
filename = ['diatC_monthly_0833-0842.nc']; ncdisp(filename)
    diatC = ncread(filename,'diatC'); %Diatom Carbon in mmol C/m^3 (size is nlon x nlat x z_t x time)   
filename = ['diazC_monthly_0833-0842.nc']; ncdisp(filename)
    diazC = ncread(filename,'diazC'); %Diazotroph Carbon in mmol C/m^3 (size is nlon x nlat x z_t x time)        
filename = ['spC_monthly_0833-0842.nc']; ncdisp(filename)
    spC = ncread(filename,'spC'); %Small Phytoplankton Carbon in mmol C/m^3 (size is nlon x nlat x z_t x time)      

%%%% Chemical tracers and gas fluxes %%%%       
filename = ['ALK_monthly_0833-0842.nc']; ncdisp(filename)
    ALK = ncread(filename,'ALK'); %Alkalinity in meq/m^3, (size is nlon x nlat x z_t x time)
filename = ['DIC_monthly_0833-0842.nc']; ncdisp(filename)
    DIC = ncread(filename,'DIC'); %DIC in mmol C/m^3 (size is nlon x nlat x z_t x time)
filename = ['DOC_monthly_0833-0842.nc']; ncdisp(filename)
    DOC = ncread(filename,'DOC'); %Dissolved Organic Carbon in mmol C/m^3 (size is nlon x nlat x z_t x time) 
filename = ['pCO2SURF_monthly_0833-0842.nc']; ncdisp(filename)
    pCO2SURF = ncread(filename,'pCO2SURF'); %Surface pCO2 in ppmv (size is nlon x nlat x time)   
    
filename = ['FG_CO2_monthly_0833-0842.nc']; ncdisp(filename)
    FG_CO2 = ncread(filename,'FG_CO2'); %Carbon Dioxide Flux in mmol DIC/m^3 cm/sec (size is nlon x nlat x time)
filename = ['FG_O2_monthly_0833-0842.nc']; ncdisp(filename)
    FG_O2 = ncread(filename,'FG_O2'); %'Oxygen Flux in mmol O2/m^3 cm/sec (size is nlon x nlat x time)
filename = ['PV_CO2_monthly_0833-0842.nc']; ncdisp(filename)
    PV_CO2 = ncread(filename,'PV_CO2'); %Piston Velocity of CO2 in cm/sec (size is nlon x nlat x time)      
filename = ['PV_O2_monthly_0833-0842.nc']; ncdisp(filename)
    PV_O2 = ncread(filename,'PV_O2'); %Piston Velocity of O2 in cm/sec (size is nlon x nlat x time)    
    
%%%% Flux of POC %%%%
filename = ['POC_FLUX_IN_monthly_0833-0842.nc']; ncdisp(filename)
    POC_FLUX_IN = ncread(filename,'POC_FLUX_IN'); %Incoming Flux of POC in mmol POC/m^3 cm/sec (size is nlon x nlat x z_t x time)   
filename = ['POC_PROD_monthly_0833-0842.nc']; ncdisp(filename)
    POC_PROD = ncread(filename,'POC_PROD'); %Production of POC in mmol POC/m^3/sec (size is nlon x nlat x z_t x time)      
filename = ['POC_REMIN_monthly_0833-0842.nc']; ncdisp(filename)
    POC_REMIN = ncread(filename,'POC_REMIN'); %Remineralization of POC in mmol POC/m^3/sec (size is nlon x nlat x z_t x time)  
    
%%%% Flux of CaCO3, SiO2, and Fe (particles and dust) %%%%    
% filename = ['bSi_form_monthly_0833-0842.nc']; ncdisp(filename)
%     bSi_form = ncread(filename,'bSi_form'); %'Production of biogenic Si', in mmol Si/m^3/sec (size is nlon x nlat x z_t x time)    
% filename = ['SiO2_FLUX_IN_monthly_0833-0842.nc']; ncdisp(filename)
%     SiO2_FLUX_IN = ncread(filename,'SiO2_FLUX_IN'); %Incoming Flux of SiO2 in mmol SiO2/m^3 cm/sec (size is nlon x nlat x z_t x time)   
% filename = ['SiO2_PROD_monthly_0833-0842.nc']; ncdisp(filename)
%     SiO2_PROD = ncread(filename,'SiO2_PROD'); %Production of SiO2 in mmol SiO2/m^3/sec (size is nlon x nlat x z_t x time)      
% filename = ['SiO2_REMIN_monthly_0833-0842.nc']; ncdisp(filename)
%     SiO2_REMIN = ncread(filename,'SiO2_REMIN'); %Remin. of SiO2 in mmol SiO2/m^3/sec (size is nlon x nlat x z_t x time)  
%     
% filename = ['CaCO3_FLUX_IN_monthly_0833-0842.nc']; ncdisp(filename)
%     CaCO3_FLUX_IN = ncread(filename,'CaCO3_FLUX_IN'); %Incoming Flux of CaCO3 in mmol CaCO3/m^3 cm/sec (size is nlon x nlat x z_t x time)   
% filename = ['CaCO3_form_monthly_0833-0842.nc']; ncdisp(filename)
%     CaCO3_form = ncread(filename,'CaCO3_form'); %Formation of CaCO3 in mmol CaCO3/m^3/sec (size is nlon x nlat x z_t x time)      
% filename = ['CaCO3_PROD_monthly_0833-0842.nc']; ncdisp(filename)
%     CaCO3_PROD = ncread(filename,'CaCO3_PROD'); %Production of Particulate CaCO3 in mmol CaCO3/m^3/sec (size is nlon x nlat x z_t x time)        
% filename = ['CaCO3_REMIN_monthly_0833-0842.nc']; ncdisp(filename)
%     CaCO3_REMIN = ncread(filename,'CaCO3_REMIN'); %Remin. of Particulate CaCO3 in mmol CaCO3/m^3/sec (size is nlon x nlat x z_t x time)     
% 
% filename = ['P_iron_FLUX_IN_monthly_0833-0842.nc']; ncdisp(filename)
%     P_iron_FLUX_IN = ncread(filename,'P_iron_FLUX_IN'); %Incoming Flux of P_iron in mmol Fe/m^3 cm/sec (size is nlon x nlat x z_t x time)   
% filename = ['P_iron_PROD_monthly_0833-0842.nc']; ncdisp(filename)
%     P_iron_PROD = ncread(filename,'P_iron_PROD'); %Production of P_iron in mmol Fe/m^3/sec (size is nlon x nlat x z_t x time)      
% filename = ['P_iron_REMIN_monthly_0833-0842.nc']; ncdisp(filename)
%     P_iron_REMIN = ncread(filename,'P_iron_REMIN'); %Remin. of P_iron in mmol Fe/m^3/sec (size is nlon x nlat x z_t x time)          
% filename = ['dust_FLUX_IN_monthly_0833-0842.nc']; ncdisp(filename)
%     dust_FLUX_IN = ncread(filename,'dust_FLUX_IN'); %Incoming Flux of dust in g/cm^2/sec (size is nlon x nlat x z_t x time)    
% filename = ['dust_REMIN_monthly_0833-0842.nc']; ncdisp(filename)
%     dust_REMIN = ncread(filename,'dust_REMIN'); %Remin. of dust in g/cm^3/sec (size is nlon x nlat x z_t x time) 
 