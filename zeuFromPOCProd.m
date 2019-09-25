function [zeu] = zeuFromPOCProd(POCprod,z_in)

% Calculates euphotic depth from POC production rates, following the
% approach of Lima et al. (2014), which defines euphotic depth as the depth
% where POC production = 1% of maximum POC production in the water column.
% Input data (intended as model output) is lat x lon x depth x time. Note
% that the function creates some anomalies when the input is raw model
% output - recommended to first interpolate in depth space using
% depth_interp_4D.

% H. I. Palevsky, January 2017

% ------------------------------------------------------------------
% INPUTS:
% z_in - depth intervals
% POCprod - must be 4D, with lat & lon as 1st and 2nd dimensions, depth as
    % 3rd dimension (z), and time as 4th dimension (could be made more
    % general in future versions)
% ------------------------------------------------------------------
% OUTPUTS:
% zeu - euphotic depth, size 1st x 2nd x 4th dimension of POCprod
% ------------------------------------------------------------------

[d1s,d2s,d3s,d4s] = size(POCprod);
zeu = NaN*ones(d1s,d2s,d4s);

for i = 1:d1s
    for j = 1:d2s
        for k = 1:d4s
            POCprodpt = squeeze(POCprod(i,j,:,k));
            if length(POCprodpt(~isnan(POCprodpt))) > 2 %only calculate if more than two non-nan values
                [C,ida1,idc] = unique(POCprodpt,'stable'); %remove rows with duplicates (which cause problems for interp1 function)
                ida2 = find(isnan(POCprodpt) == 0); %remove rows with nans
                ida = intersect(ida1,ida2,'stable'); %index of rows without nans or duplicates
                if length(ida) > 2 %check if > 2 points remain for interpolation
                    zeu(i,j,k) = interp1(POCprodpt(ida),z_in(ida),0.01*max(POCprodpt));
                end
            end
        end
    end
end

end