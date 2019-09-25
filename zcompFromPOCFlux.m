function [zcomp, POCflux_zcomp] = zcompFromPOCFlux(z_in, POCflux_in)

% Calculates particle compensation depth from POC flux rates, where the
% particle compensation depth is defined as the depth of maximum POC flux
% in the water column.
% Input data (intended as model output) is lat x lon x depth x time. Note
% that the maximum POC flux is not altered (can only be slightly reduced)
% by interpolating to a finer depth grid, so this yields optimal results
% when using original model output.

% H. I. Palevsky, January 2017

% ------------------------------------------------------------------
% INPUTS:
% z_in - depth intervals
% POCflux_in - must be 4D, with lat & lon as 1st and 2nd dimensions, depth as
    % 3rd dimension (z), and time as 4th dimension (could be made more
    % general in future versions)
% ------------------------------------------------------------------
% OUTPUTS:
% zcomp - particle compensation depth, size 1st x 2nd x 4th dimension of
    % POCflux_in
% POCflux_zcomp - POC flux rate at zcomp
% ------------------------------------------------------------------

[d1s,d2s,d3s,d4s] = size(POCflux_in);
zcomp = NaN*ones(d1s,d2s,d4s);

[POCflux_zcomp, idzcomp] = max(POCflux_in,[],3);
    POCflux_zcomp = squeeze(POCflux_zcomp);
    idzcomp = squeeze(idzcomp);
 
 for i = 1:d1s
     for j = 1:d2s
         for k = 1:d4s
             zcomp(i,j,k) = z_in(idzcomp(i,j,k));
         end
     end
 end
 
end