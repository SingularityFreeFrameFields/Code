function [narrowBandNew] = repairMask(narrowBand)
 nonZerosAroundMe = conv2(narrowBand,[0 1 1; 1 0 1; 1 1 0]);
 nonZerosAroundMe = nonZerosAroundMe(2:end-1,2:end-1);
 narrowBandNew = narrowBand;
 narrowBandNew(nonZerosAroundMe>=5) = 1;
end