function [inds,gridpt] = findNearestPoint(TLAT,TLONG,stnlat,stnlon);

%Function to find the nearest grid point to an input point on an irregular
%grid (i.e. CSM grid)
latdif = abs(TLAT - stnlat);
londif = abs(TLONG - stnlon);
diftot = latdif + londif;
[Y,I] = min(diftot);
[Y2,I2] = min(Y);
inds = [I(I2),I2];
gridpt = double([TLAT(inds(1),inds(2)),TLONG(inds(1),inds(2))]);

end
