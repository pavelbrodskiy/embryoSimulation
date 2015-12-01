function [newBoundaries, engrailed] = proliferationSimulation(cellBoundaries, divisionOccurs, engrailed)

boundID = find(cellBoundaries);
dividingCells = divisionOccurs;
newBoundaries = cellBoundaries;

if ~isempty(dividingCells)
    for i = length(dividingCells):-1:1
        newBoundary = uint16(mean([boundID(dividingCells(i)),boundID(dividingCells(i)+1)]));
        newBoundaries(newBoundary) = 1;    
        engrailed = [engrailed(1:i) engrailed(i) engrailed((i+1):end)];
    end
end

end

