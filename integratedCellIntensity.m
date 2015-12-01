function [ cellIntensities, cellPlot, cellSize ] = integratedCellIntensity( cellBounds, intensity )

boundID = find(cellBounds);

cellNumber = length(boundID)-1;

for i = cellNumber:-1:1
    cellIntensities(i) = sum(intensity(boundID(i):boundID(i+1)))/(boundID(i+1)-boundID(i));
    cellSize(i) = boundID(i+1)-boundID(i);
    cellPlot(boundID(i):boundID(i+1))=cellIntensities(i);
end

end

