function [data_weightedAnnualMean] = weightedAnnualMean(dataForMean, dataForWeighting, year)

% -----------------------------------------------------
% INPUTS:
% dataForMean: size lat x lon x year
% dataForWeighting: size lat x lon x year
% year: identifier for matching when dataForMean includes multiple years
% ----------------------------------------------------
% OUTPUTS:
% data_weightedAnnualMean: size lat x lon x number of unique values in year
% ----------------------------------------------------

[d1s, d2s, numyrs] = size(dataForMean);
[m, n, p] = size(dataForWeighting);
if numyrs ~= length(year)
   error('3rd dimension of dataForMean must be same length as year')
end
if (m~=d1s | n~=d2s | p~=numyrs)
    error('dataForMean and dataForWeighting must be same size')
end

yrslist = unique(year);
data_weightedAnnualMean = NaN*(ones(d1s,d2s,length(yrslist)));


for j = 1:d1s
    for k = 1:d2s
        for i = 1:length(yrslist)
            idyr = find(year == yrslist(i));
                weight = squeeze(dataForWeighting(j,k,idyr))./(length(idyr)*mean(dataForWeighting(j,k,idyr),3));
            data_weightedAnnualMean(j,k,i) = sum(squeeze(dataForMean(j,k,idyr)).*weight);
        end
    end
end

end