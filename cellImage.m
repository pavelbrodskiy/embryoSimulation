function [ image ] = cellImage( cellBounds, engrailed )

cellCoords = find(cellBounds);
cellNumber = length(cellCoords) - 1;

for i = cellNumber:-1:1
    image(1:300, cellCoords(i):cellCoords(i+1), 2) = engrailed(i);
end

image(:,find(cellBounds),:) = 0;
image(:,:,3)=0;

image = double(image);

end

