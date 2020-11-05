function data = normRunTriplet(data)
for i = 1:numel(data.ori1)
    %average runs in cartesian space
    dataMean = mean(cat(3,data.ori1{i},data.ori2{i},data.plaid{i}),3);
    %normalize runs in cartesian space
    data.ori1{i} = data.ori1{i} - dataMean;
    data.ori2{i} = data.ori2{i} - dataMean;
    data.plaid{i} = data.plaid{i} - dataMean;
end
