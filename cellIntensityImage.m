function [ output_args ] = cellIntensityImage( cellBounds, intensity )
boundID = find(cellBounds);

cellNumber = length(boundID) - 1;

for i = cellNumber:-1:1
    cellIntensities(i) = sum(intensity(boundID(i):boundID(i+1)));
end

end


