function [newBoundaries, engrailed] = deathSimulation(cellBoundaries, deathOccurs, engrailed)

boundID = find(cellBoundaries);
newBoundaries = cellBoundaries;

if ~isempty(deathOccurs)
    for i = length(deathOccurs):-1:1
        newBoundary = uint16(mean([boundID(deathOccurs(i)),boundID(deathOccurs(i)+1)]));
        newBoundaries(newBoundary) = 1;
        newBoundaries(boundID(deathOccurs(i))) = 0;
        newBoundaries(boundID(deathOccurs(i)+1)) = 0;
        engrailed(i) = [];
    end
end

end

