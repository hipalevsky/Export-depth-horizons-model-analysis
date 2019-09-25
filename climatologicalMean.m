function [data_clim] = climatologicalMean(data)

% Calculates climatological mean values of 3 or 4 dimensional data, where
% the final dimension is time and has 12 monthly entries in each year

% H. I. Palevsky, January 2017

% ------------------------------------------------------------------
% INPUTS:
% data - size can be 3 or 4 dimensions, with time in last dimension and is
    % structured as 12 repeating months (so if last dimension has length 12, it
    % includes 1 year of data; if it has length 120, it has 10 years of data)
% ------------------------------------------------------------------
% OUTPUTS:
% data_clim - climatological mean monthly values of data
% ------------------------------------------------------------------

N = ndims(data);
numtimes = size(data,N); %time must be in last dimension

for i = 1:12
    idmon = [i:12:numtimes];
    if N == 3
        data_clim(:,:,i) = mean(data(:,:,idmon),3);
    elseif N == 4
        data_clim(:,:,:,i) = mean(data(:,:,:,idmon),4);
    end
end

end