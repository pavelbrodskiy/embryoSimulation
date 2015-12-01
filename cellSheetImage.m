function [ cellImage ] = cellSheetImage( boundaries, signal )

cellBounds = 100 * (find(boundaries)-1)+1;
cellNumber = length(cellBounds) - 1;

for i = cellNumber:-1:1
    midpoints(i) = mean(cellBounds(i),cellBounds(i+1));
    width(i) = cellBounds(i+1) - cellBounds(i);
    cellImage(1,cellBounds(i):cellBounds(i+1)) = signal(midpoints(i));
end

width(width > 1000) = 1000;

cellImage = repmat(cellImage,[5000,1]);
cellImage(:,:,2) = 0;
cellImage(:,:,3) = 0;

end

