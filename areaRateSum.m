function [sum_glob, sum_reg, area_reg] = areaRateSum(rate_grid,latdiv,TLAT,TAREA,REGION_MASK)
%Calculate the total rate of the parameter in rate_grid (for each grid
%cell), both globally and for each zonal region defined by bands with
%minimim and maximum latitudes given in latdiv
%--------------------------------------------------------
%INPUT:
%rate_grid: globally-gridded rates (i.e. on CSM grid) in mol C m-2 yr-1
%latdiv: minimum and maximum latitude bands for zonal region divisions
%TLAT: latitude for each grid cell (from CSM)
%TAREA: grid cell areas from CSM (in cm^2)
%REGION_MASK: mask to identify ocean cells to use from CSM
%
%OUTPUT:
%sum_glob = global rate in Pg C yr-1
%sum_reg = rate for each zonal region in Pg C yr-1
%area_reg = area for each region in m^2


    cmToMeters = 1/100;
    molToPg = 12*10^-15;
numreg = length(latdiv) - 1;
A = rate_grid.*TAREA*(cmToMeters)^2*molToPg;
sum_glob = nansum(nansum(A));
for i = 1:numreg
    reg_id = find(TLAT > latdiv(i) & TLAT <= latdiv(i+1) & REGION_MASK > 0);
    sum_reg(i) = nansum(A(reg_id));
    area_reg(i) = nansum(TAREA(reg_id))*(cmToMeters)^2;
end

end