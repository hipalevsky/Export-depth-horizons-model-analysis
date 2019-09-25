function [data_new] = depth_interp_4D(z_orig, z_new, data)

% Passes in 4D data (i.e. lon x lat x depth x time) with depth in 3rd
% dimension and interpolates from the original depth grid z_orig to a new
% depth grid z_new. Currently written only using linear interpolation and
% where depth must be in the 3rd dimension.
%
% H. I. Palevsky, January 2017

% ------------------------------------------------------------------
% INPUTS:
% z_orig - original depth intervals
% z_new - new depth grid to interpolate onto
% data - must be 4D, with depth as 3rd dimension (could be made more
    % general in future versions)
% ------------------------------------------------------------------
% OUTPUTS:
% data_new - 1st, 2nd, and 4th dimensions remain same size as before, now
% interpolated in 3rd dimension from from z_orig to z_new
% ------------------------------------------------------------------

[d1s,d2s,d3s,d4s] = size(data);
data_new = NaN*ones(d1s,d2s,length(z_new),d4s); %initialize array to hold output

for i = 1:d4s
    A = interp1(z_orig,shiftdim(data(:,:,:,i),2),z_new);
    data_new(:,:,:,i) = shiftdim(A,1);
end

end