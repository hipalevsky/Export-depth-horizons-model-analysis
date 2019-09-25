function [vals] = depthCriteriaVals(data, z_in, criteria)

% Calculates values of model data (size lat x lon x depth x time)
% interpolated to depth criteria (input in cell array criteria, where each
% criterion is size lat x lon x time).

% H. I. Palevsky, January 2017

% ------------------------------------------------------------------
% INPUTS:
% data - size lat x lon x depth x time
% z_in - depth intervals
% criteria - cell array of depth criteria (same units as z_in) to use to
    % interpolate data. Can include any number of cells (i.e. criteria),
    % but each criterion must be size lat x lon x time.
% ------------------------------------------------------------------
% OUTPUTS:
% vals - values of data, interpolated onto the depths specified in
    % criteria, with the same size and structure as criteria.
% ------------------------------------------------------------------

[d1s,d2s,d3s,d4s] = size(data);
numcriteria = length(criteria);

%Initialize output
for i = 1:numcriteria
    vals{i} = NaN*ones(d1s,d2s,d4s);
end

for i = 1:d1s
    for j = 1:d2s
        for k = 1:d4s
            for c = 1:numcriteria
                critpt(c) = criteria{c}(i,j,k); %Pull out all depth criteria for given time-location point.
            end
            interppt = interp1(z_in,squeeze(data(i,j,:,k)),critpt); %Calculate data values for all depth criteria.
            for c = 1:numcriteria
                vals{c}(i,j,k) = interppt(c); %Assign interpolated values to output for the given location-time point for the cell corresponding to the given depth criterion.
            end
        end
    end
end


end