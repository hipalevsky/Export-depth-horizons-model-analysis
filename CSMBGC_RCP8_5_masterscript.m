%% Wrapper script for all components of CSMBGC RCP8.5 analysis

%% Run startup script to bring in relevant toolboxes
% startup

%% Add path for needed output
addpath('C:/Users/palevsky/Dropbox/MATLAB/CSMBGC data processing/RCP 8_5 processed output')

%% Extraction of annual data and calculation of climatological means over 20 year periods
% Follows the model in CSMBGC_exportanalysis_RCP8_5 (for varying year
% ranges set in line 13)

%% Comparison of beginning and end of the 21st century (20 yr climatological means)

% 1) Load previously calculated climatological mean data from beginning and
% end of century and prepare to use for Taylor decomposition calculations
% and comparisons
% CSMBGC_RCP8_5_changeOverTime

% 2) Taylor decomposition calculations
% CSMBGC_RCP8_5_TaylorDecomposition

% 3) Regrid previous output (from last two scripts) and save for plotting
% CSMBGC_RCP8_5_changeOverTime_regrid

% 4) Script to plot global maps of climatological data compared between the
%beginning and end of the 21st century
CSMBGC_RCP8_5_changeOverTime_plotting

%% Time-series site analysis
% Script used to extract and process data:
% CSMBGC_RCP8_5_TimeSeriesSites

% Script to interpret and plot the time series site data
CSMBGC_RCP8_5_TimeSeriesSites_plotting

%% Comparison of MLDmax with observations from WOA
MLDmax_model_WOA_compare