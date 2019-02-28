function quantid = checkQuant(outsiz,inputloc)
bins = outsiz/2 < inputloc;
quantid = sum(bins.*(2.^[0 1 2]))+1;


